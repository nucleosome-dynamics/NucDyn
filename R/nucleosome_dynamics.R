#!/usr/bin/env Rscript

setMethod("nucleosomeDynamics", signature(setA="IRanges", setB="IRanges"),
    function(setA, setB, equalSize, ...) {
        sets <- list(setA, setB)
        myDyn <- .nucleosomeDynamics(sets, equalSize=equalSize, ...)
        # we wrap it in a list so that buildNucDyn behaves as expected
        myDyn <- .buildNucDyn(list("*"=myDyn), equalSize)
        myDyn
    }
)

setMethod("nucleosomeDynamics", signature(setA="RangedData", setB="RangedData"),
    function(setA, setB, equalSize, mc.cores=1, ...) {
        sets <- list(setA, setB)

        # do it for every chromosome separately
        chrs <- unique(unlist(lapply(sets, function(set) levels(space(set)))))
        myDyn <- .xlapply(
            chrs,
            function(chr) {
                message(paste("Starting", chr))
                dyn <- .nucleosomeDynamics(lapply(sets, function(x) x[chr]),
                                           equalSize=equalSize, ...)
                message(paste(chr, "done"))
                dyn
            },
            mc.cores=mc.cores
        )
        names(myDyn) <- chrs

        # join it back together into a single list of pairs of GRanges
        myDyn <- .buildNucDyn(myDyn, equalSize)

        myDyn
    }
)

setMethod("nucleosomeDynamics", signature(setA="GRanges", setB="GRanges"),
    function(setA, setB, equalSize, mc.cores=1, ...) {
        sets <- list(setA, setB)

        splitted <- lapply(
            sets,
            function(x)
                GenomicRanges::split(x, seqnames(x))
        )
        chrs <- unique(unlist(lapply(splitted, names)))

        # do it for every chromosome separately
        myDyn <- .xlapply(
            chrs,
            function(chr) {
                message(paste("Starting", chr))

                dyn <- NucDyn:::.nucleosomeDynamics(
                    lapply(splitted,
                           function(x) ranges(x[[chr]])),
                    equalSize=equalSize,
                    ...
                )

                message(paste(chr, "done"))
                dyn
            },
            mc.cores=mc.cores
        )
        names(myDyn) <- chrs

        # join it back together into a single list of pairs of GRanges
        myDyn <- .buildNucDyn(myDyn, equalSize)

        myDyn
    }
)

.buildNucDyn <- function(dyn, equalSize)
{  # Build a NucleosomeDynamics object
    message("Combining the calculations from different chromosomes...")
    NucDyn(set.a=.getGRLs(dyn, 1, equalSize),
           set.b=.getGRLs(dyn, 2, equalSize))
}

.getGRLs <- function(dyn, i, equalSize)
{  # Build a GRangesList for a given set

    if (equalSize) {
        readTypes <- c("originals", "coinciding", "left.shifts",
                       "right.shifts", "indels", "unpaired")
    } else {
        readTypes <- c("originals", "coinciding", "same.start",
                       "same.end", "containedA", "containedB",
                       "left.shifts", "right.shifts", "indels",
                       "unpaired")
    }

    setLs <- lapply(dyn, function(x) lapply(x, function(y) y[[i]]))

    lens <- lapply(setLs, function(x) sapply(x, length))

    gr <- GRanges(
        seqnames = Rle(names(setLs), sapply(lens, sum)),
        ranges = do.call("c", unname(unlist(setLs))),
        type = Rle(sapply(setLs, names), unlist(lens))
    )
    grLs <- GenomicRanges::split(gr, gr$type)  # split by type

    for (type in readTypes) {  # add possibly missing types
        if (!type %in% names(grLs)) {
            grLs[[type]] <- GRanges()
        }
    }

    grLs[readTypes]  # keep them in the wanter order
}

.nucleosomeDynamics <- function(mySets, maxLen=170, roundPow=5, equalSize,
                                readSize=140, maxDiff=74)
{
    mySets <- lapply(mySets, .toIRanges)
    originals <- mySets

    # remove reads that are too long
    mySets <- lapply(mySets, function(x) .rmLongReads(x, maxLen=maxLen))
    # make the mean read length of both sets equal
    if (equalSize) {
        #mySets <- .stdSetLen(mySets)
        ## remove reads that are too long (again)
        #mySets <- lapply(mySets, function(x) .rmLongReads(x, maxLen=maxLen))

        mySets <- lapply(mySets, function(x) .setSizeTo(x, readSize=readSize))
        mySets <- lapply(mySets, IRanges::sort)  # keep them sorted

        # reads considered to be the same with a small variance distance allowed
        subsetList <- .equalsAtDist(mySets, maxDist=5)
        newSets <- .separateGroups(mySets, subsetList)
        equalReads <- newSets$matches
        mySets <- newSets$rest

    } else {
        # round the reads to powers of 5
        mySets <- lapply(mySets, function(x) .setRounder(x, pow=roundPow))
        mySets <- lapply(mySets, IRanges::sort)  # keep them sorted

        # pairs that start and end at the same position
        subsetList <- .equals(mySets)
        newSets <- .separateGroups(mySets, subsetList)
        equalReads <- newSets$matches
        mySets <- newSets$rest

        # pairs that start at the same position but end at a different one
        subsetList <- .sameStart(mySets)
        newSets <- .separateGroups(mySets, subsetList)
        sameStartReads <- newSets$matches
        mySets <- newSets$rest

        # pairs that end at the same position but start at a different one
        subsetList <- .sameEnd(mySets)
        newSets <- .separateGroups(mySets, subsetList)
        sameEndReads <- newSets$matches
        mySets <- newSets$rest

        subsetListA <- .contained(mySets)
        subsetListB <- rev(.contained(rev(mySets)))

        containedReadsA <- .separateGroups(mySets, subsetListA)$matches
        containedReadsB <- .separateGroups(mySets, subsetListB)$matches

        containedSubset <- mapply(
            function(x, y)
                as.numeric(as.logical(x) | as.logical(y)),
            subsetListA,
            subsetListB
        )

        mySets <- .separateGroups(mySets, containedSubset)$rest
    }

    # preliminar shift search
    preShiftSubsets <- .preShifts(mySets, maxDiff)
    newSets <- .separateShifts(mySets, preShiftSubsets)
    preLeftShiftReads <- newSets$left
    preRightShiftReads <- newSets$right

    shiftClust <- .clusterizeShifts(preLeftShiftReads[[1]], preRightShiftReads[[1]])

    shiftSubsets <- .shifts(mySets, shiftClust, maxDiff)
    newSets <- .separateShifts(mySets, shiftSubsets)
    leftShiftReads <- newSets$left
    rightShiftReads <- newSets$right
    mySets <- newSets$rest

    equalNumbered <- .normSize(mySets)
    discarded <- equalNumbered$removed
    mySets <- equalNumbered$rest

    if (equalSize) {
        returnLs <- list(originals    = originals,
                         coinciding   = equalReads,
                         left.shifts  = leftShiftReads,
                         right.shifts = rightShiftReads,
                         indels       = mySets,
                         unpaired     = discarded)
    } else {
        returnLs <- list(originals    = originals,
                         coinciding   = equalReads,
                         same.start   = sameStartReads,
                         same.end     = sameEndReads,
                         containedA   = containedReadsA,
                         containedB   = containedReadsB,
                         left.shifts  = leftShiftReads,
                         right.shifts = rightShiftReads,
                         indels       = mySets,
                         unpaired     = discarded)
    }

    return(returnLs)
}
