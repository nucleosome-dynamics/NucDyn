#!/usr/bin/env Rscript

setMethod(
    "nucleosomeDynamics",
    signature(setA="IRanges", setB="IRanges"),
    function(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140,
             maxDiff=74, minDiff=10) {
        sets <- list(setA, setB)
        myDyn <- .nucleosomeDynamics(mySets=sets,
                                     maxLen=maxLen,
                                     roundPow=roundPow,
                                     equalSize=equalSize,
                                     readSize=readSize,
                                     maxDiff=maxDiff,
                                     minDiff=minDiff)
        # we wrap it in a list so that buildNucDyn behaves as expected
        myDyn <- .buildNucDyn(list("*"=myDyn), equalSize)
        myDyn
    }
)

setMethod(
    "nucleosomeDynamics",
    signature(setA="RangedData", setB="RangedData"),
    function(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140,
             maxDiff=74, minDiff=10, mc.cores=1) {
        sets <- list(setA, setB)

        # do it for every chromosome separately
        chrs <- unique(unlist(lapply(sets, function(set) levels(space(set)))))
        myDyn <- .xlapply(
            chrs,
            function(chr) {
                message(paste("Starting", chr))
                dyn <- .nucleosomeDynamics(mySets=lapply(sets, "[", chr),
                                           maxLen=maxLen,
                                           roundPow=roundPow,
                                           equalSize=equalSize,
                                           readSize=readSize,
                                           maxDiff=maxDiff,
                                           minDiff=minDiff)
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

setMethod(
    "nucleosomeDynamics",
    signature(setA="GRanges", setB="GRanges"),
    function(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140,
             maxDiff=74, minDiff=10, mc.cores=1) {
        sets <- list(setA, setB)

        splitted <- lapply(
            sets,
            function(x)
                split(x, seqnames(x))
        )
        chrs <- unique(unlist(lapply(splitted, names)))

        # do it for every chromosome separately
        myDyn <- .xlapply(
            chrs,
            function(chr) {
                message(paste("Starting", chr))

                dyn <- .nucleosomeDynamics(
                    mySets=lapply(splitted, function(x) ranges(x[[chr]])),
                    maxLen=maxLen,
                    roundPow=roundPow,
                    equalSize=equalSize,
                    readSize=readSize,
                    maxDiff=maxDiff,
                    minDiff=minDiff
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
        readTypes <- c("originals", "coinciding", "small.shifts",
                       "left.shifts", "right.shifts", "indels",
                       "unpaired")
    } else {
        readTypes <- c("originals", "coinciding", "same.start", "same.end",
                       "containedA", "containedB", "small.shifts",
                       "left.shifts", "right.shifts", "indels", "unpaired")
    }

    setLs <- lapply(dyn, function(x) lapply(x, function(y) y[[i]]))

    lens <- lapply(setLs, function(x) sapply(x, length))

    gr <- GRanges(
        seqnames = Rle(names(setLs), sapply(lens, sum)),
        ranges = do.call("c", unname(unlist(setLs))),
        type = Rle(sapply(setLs, names), unlist(lens))
    )
    grLs <- split(gr, gr$type)  # split by type

    for (type in readTypes) {  # add possibly missing types
        if (!type %in% names(grLs)) {
            grLs[[type]] <- GRanges()
        }
    }

    grLs[readTypes]  # keep them in the wanted order
}

.nucleosomeDynamics <- function(mySets, maxLen=170, roundPow=5, equalSize,
                                readSize=140, maxDiff=74, minDiff=10)
{
    mySets <- lapply(mySets, .toIRanges)
    #originals <- mySets

    # remove reads that are too long
    mySets <- lapply(mySets, .rmLongReads, maxLen=maxLen)
    # make the mean read length of both sets equal
    if (equalSize) {
        #mySets <- .stdSetLen(mySets)
        ## remove reads that are too long (again)
        #mySets <- lapply(mySets, function(x) .rmLongReads(x, maxLen=maxLen))

        mySets <- lapply(mySets, .setSizeTo, readSize=readSize)
        #mySets <- lapply(mySets, IRanges::unique)
        mySets <- lapply(mySets, sort)  # keep them sorted

        originals <- mySets

        # reads considered to be the same with a small variance distance allowed
        subsetList <- .equalsAtDist(mySets, maxDist=5)
        newSets <- .separateGroups(mySets, subsetList)
        equalReads <- newSets$matches
        mySets <- newSets$rest

    } else {
        # round the reads to powers of 5
        mySets <- lapply(mySets, .setRounder, pow=roundPow)
        #mySets <- lapply(mySets, IRanges::unique)
        mySets <- lapply(mySets, sort)  # keep them sorted

        originals <- mySets

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

        mySets <- mapply(
            function(set, subA, subB)
                set[!(as.logical(subA) | as.logical(subB))],
            mySets,
            subsetListA,
            subsetListB
        )

    }

    newSets <- .shifts(mySets, max.dist=maxDiff, min.dist=0)

    postLeft <- .applyDistThresh(newSets$left, minDiff)
    postRight <- .applyDistThresh(newSets$right, minDiff)

    smallShiftReads <- mapply(c,
                              postLeft$smallShifts,
                              postRight$smallShifts,
                              SIMPLIFY=FALSE)

    leftShiftReads <- postLeft$properShifts
    rightShiftReads <- postRight$properShifts

    mySets <- newSets$rest

    #equalNumbered <- .normSize(mySets)
    #discarded <- equalNumbered$removed
    #mySets <- equalNumbered$rest
    discarded <- list(IRanges(), IRanges())

    if (equalSize) {
        returnLs <- list(originals    = originals,
                         coinciding   = equalReads,
                         small.shifts = smallShiftReads,
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
                         small.shifts = smallShiftReads,
                         left.shifts  = leftShiftReads,
                         right.shifts = rightShiftReads,
                         indels       = mySets,
                         unpaired     = discarded)
    }

    return(returnLs)
}
