#' Run nucleosomeDynamics to compare two sets.
#'
#' This is the main function in nucleosomeDynamics. It allows to compare the
#' reads of two NGS experiments of nucleosome coverage.
#'
#' The aim of NucleosomeDynamics is to infer "movement" (with direction and
#' magnitude) of the reads between two reference nucleosome maps. In contrast
#' with a simple coverage difference, NucleosomeDynamics can tell how the reads
#' change between two different experiments. This is useful to analyze regions
#' where fine regulatory role of the nucleosomes is suspected to happen.
#'
#' This method is based on the idea that reads in a reference state (ref1)
#' should match those in another reference state (ref2) after applying a few
#' shifts and/or indels. Both ref1 and ref2 need to be experimental nucleosome
#' maps, either from the same sample with different conditions or from
#' different samples.
#'
#' Then, we look for a match to each read in ref1 in another read of ref2 using
#' a specific deffinition for what a "match" is. After all possible matches
#' have been found, we set those reads apart and we look for matches in the
#' remaining reads using a different definition for what a "match" is.  Then,
#' to account for the possibility that one sample has a higher coverage than
#' the other, randomly picked reads are removed from the dataset with more
#' reads to that they are the same size.  After this, all the remaining reads
#' are considered indels.
#'
#' Since they are tried sequentially, the definition of each type of match
#' implies that the previous definitions tried are do not hold.  The different
#' types of matches, in the order in which they are tried are:
#'
#' * Coinciding: Reads that start and end in the exact same position in both
#'   sets.
#' * Same start: Reads that start in the same position in both sets but end at
#'   a different one.
#' * Same end: Reads that end in the same position in both sets but start at a
#'   different one.
#' * Contained: Reads from one set that are contained or contain reads from the
#'   other set. For a read to be contained by another, it has to start at a
#'   more upstream position but end in a more downstream position than the
#'   second read.
#' * Shifts: Reads whose dyads are at a maximum distance of `maxDist` (default
#'   is 74 bp.).
#'
#' If `equalSize=TRUE`, all reads are forced to the same size, to the match
#' types "Same start", "Same end" and "Contained" do not apply.
#'
#' In an attempt to find the optimum pairing for the shifts, we use dynamic
#' programming approach. It is done in such a way that the score is inversely
#' proportional to the dyad distance (to favor shortest possible distance), but
#' with a very high gap penalty (to favor more pairs at longer distance rather
#' than less closer pairs) and a penalty of `-Inf` for distances higher than
#' `maxDist` so that those don't happen at all.
#'
#' @param setA Reads of the first experiment in `IRanges`, `RangedData`, or
#'   `GRanges`. The format should the same of the one in setB.
#' @param setB Reads of the second experiment in `IRanges`, `RangedData`, or
#'   `GRanges. The format should the the same of the one in setA.
#' @param maxLen Reads longer than this number will be filtered out.
#' @param equalSize If set to `TRUE`, all sets will be set to the same length,
#'   conserving their original dyad position. Use it if the reads in your sets
#'   have differences in length (ie, due to differences in the digestion) that
#'   you are not interested in.
#' @param readSize Length to which all reads will be set in case `equalSize` is
#'   `TRUE`. It is ignored when `equalSize` is set to `FALSE`.
#' @param maxDiff Maximum distance between the dyads of two reads that allows
#'   them to still be considered a "shift".
#' @param minDiff Minimum distance between the dyads of two reads that allows
#'   them to still be considered a "shift".
#' @param mc.cores If `parallel` support, the number of cores available. This
#'   option is only used if the provided sets are from more than one
#'   chromosome.
#'
#' @return An object of class [NucDyn-class].
#'
#' @examples
#'     \donttest{data(readsG2_chrII)}
#'     \donttest{data(readsM_chrII)}
#'     \donttest{dyn <- nucleosomeDynamics(setA=readsG2_chrII, setB=readsM_chrII)}
#'
#' @author Ricard Illa,
#'         Diana Buitrago \email{diana.buitrago@@irbbarcelona.org}
#' @keywords manip
#' @rdname nucleosomeDynamics
#' @export nucleosomeDynamics
#'
setGeneric(
    "nucleosomeDynamics",
    function(setA, setB, maxLen=170, equalSize=FALSE, readSize=140,
             maxDiff=74, minDiff=10, mc.cores=1)
        standardGeneric("nucleosomeDynamics")
)

#' @rdname nucleosomeDynamics
setMethod(
    "nucleosomeDynamics",
    signature(setA="IRanges", setB="IRanges"),
    function(setA, setB, maxLen=170, equalSize=FALSE, readSize=140,
             maxDiff=74, minDiff=10) {
        sets <- list(setA, setB)
        myDyn <- .nucleosomeDynamics(mySets=sets,
                                     maxLen=maxLen,
                                     equalSize=equalSize,
                                     readSize=readSize,
                                     maxDiff=maxDiff,
                                     minDiff=minDiff)
        # we wrap it in a list so that buildNucDyn behaves as expected
        myDyn <- .buildNucDyn(list("*"=myDyn), equalSize)
        myDyn
    }
)

#' @rdname nucleosomeDynamics
#' @importMethodsFrom IRanges space
setMethod(
    "nucleosomeDynamics",
    signature(setA="RangedData", setB="RangedData"),
    function(setA, setB, maxLen=170, equalSize=FALSE, readSize=140,
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

#' @rdname nucleosomeDynamics
#' @importMethodsFrom GenomeInfoDb seqnames
#' @importMethodsFrom IRanges ranges
setMethod(
    "nucleosomeDynamics",
    signature(setA="GRanges", setB="GRanges"),
    function(setA, setB, maxLen=170, equalSize=FALSE, readSize=140,
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

#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom S4Vectors Rle
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

    setLs <- lapply(dyn, lapply, `[[`, i)
    lens <- lapply(setLs, sapply, length)

    gr <- GRanges(
        seqnames = Rle(names(setLs), sapply(lens, sum)),
        ranges   = do.call("c", unname(unlist(setLs))),
        type     = Rle(sapply(setLs, names), unlist(lens))
    )
    grLs <- S4Vectors::split(gr, gr$type)  # split by type

    for (type in readTypes) {  # add possibly missing types
        if (!type %in% names(grLs)) {
            grLs[[type]] <- GRanges()
        }
    }

    grLs[readTypes]  # keep them in the wanted order
}

#' @importFrom IRanges IRanges
.nucleosomeDynamics <- function(mySets, maxLen=170, equalSize,
                                readSize=140, maxDiff=74, minDiff=10)
{
    mySets <- lapply(mySets, .toIRanges)
    #originals <- mySets

    # remove reads that are too long
    mySets <- lapply(mySets, .rmLongReads, maxLen=maxLen)
    # make the mean read length of both sets equal
    if (equalSize) {
        mySets <- lapply(mySets, .setSizeTo, readSize=readSize)
        mySets <- lapply(mySets, sort)  # keep them sorted

        originals <- mySets

        # reads considered to be the same with a small variance distance
        # allowed
        subsetList <- equals_at_dist(mySets, max_dist=5)
        newSets <- .separateGroups(mySets, subsetList)
        equalReads <- newSets$matches
        mySets <- newSets$rest

    } else {
        mySets <- lapply(mySets, sort)  # keep them sorted

        originals <- mySets

        # pairs that start and end at the same position
        subsetList <- equals(mySets[[1]], mySets[[2]])
        newSets <- .separateGroups(mySets, subsetList)
        equalReads <- newSets$matches
        mySets <- newSets$rest

        # pairs that start at the same position but end at a different one
        subsetList <- same_start(mySets[[1]], mySets[[2]])
        newSets <- .separateGroups(mySets, subsetList)
        sameStartReads <- newSets$matches
        mySets <- newSets$rest

        # pairs that end at the same position but start at a different one
        subsetList <- same_end(mySets[[1]], mySets[[2]])
        newSets <- .separateGroups(mySets, subsetList)
        sameEndReads <- newSets$matches
        mySets <- newSets$rest

        subsetListA <- contained(mySets[[1]], mySets[[2]])
        subsetListB <- rev(contained(mySets[[2]], mySets[[1]]))

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

    newSets <- shifts(mySets[[1]], mySets[[2]], max.dist=maxDiff, min.dist=0)

    postLeft <- .applyDistThresh(newSets$left, minDiff)
    postRight <- .applyDistThresh(newSets$right, minDiff)

    smallShiftReads <- mapply(c,
                              postLeft$smallShifts,
                              postRight$smallShifts,
                              SIMPLIFY=FALSE)

    leftShiftReads <- postLeft$properShifts
    rightShiftReads <- postRight$properShifts

    mySets <- newSets$rest

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
