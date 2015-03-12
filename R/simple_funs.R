# wrapper to choose between lapply and mclapply accordingly
.xlapply <- function(X, FUN, mc.cores=1, ...)
{
    if (mc.cores > 1) {
        succ.mc <- "parallel" %in% loadedNamespaces()
        if (!succ.mc) {
            warning("'parallel' library not available, switching to mc.cores=1")
            return(lapply(X=X, FUN=FUN, ...))
        } else {
            return(mclapply(X=X, FUN=FUN, mc.cores=mc.cores, ...))
        }
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}

# round a number to a power of 5 to make things more discrete
.roundToPower <- function(n, pow=5)
    round(n/pow) * pow

.setRounder <- function(set, pow=5)
{  # make a set more discrete
    start(set) <- .roundToPower(start(set), pow=pow)
    end(set) <- .roundToPower(end(set), pow=pow)
    set
}

.rmLongReads <- function(set, maxLen=170)
{  # remove reads longer than a specified maximum length
    set[width(set) < maxLen, ]
}

.stdSetLen <- function(sets)
{  # make the median read length of both sets roughly equal, by enlarging
   # the reads in the set with the smallest median length
    set1 <- sets[[1]]
    set2 <- sets[[2]]
    lenDiff <- median(width(set1)) - median(width(set2))
    halfLen <- round(abs(lenDiff)/2)

    if (lenDiff > 0) {  # median read length is higher in set1
        start(set2) <- start(set2) - halfLen
        end(set2) <- end(set2) + halfLen
    } else if (lenDiff < 0) {  # higher in set2
        start(set1) <- start(set1) - halfLen
        end(set1) <- end(set1) + halfLen
    }

    list(set1, set2)
}

.normSize <- function (sets)
{
    sizeA <- length(sets[[1]])
    sizeB <- length(sets[[2]])

    sDiff <- sizeA - sizeB

    if (sDiff > 0) {  # first one is bigger
        biggest <- sets[[1]]
        bigSize <- sizeA
    } else if (sDiff < 0) {  # the second one is bigger
        biggest <- sets[[2]]
        bigSize <- sizeB
    } else {  # the improbable case where they are exactly the same size
        return(list(removed=list(NULL, NULL), rest=sets))
    }

    toRm <- sample(1:bigSize, abs(sDiff))
    normalized <- biggest[-toRm, ]
    removed <- biggest[toRm, ]
    # Keep them sorted. Blame my OCD.
    sortIdxs <- sort(start(removed), index.return=TRUE)[[2]]
    removed <- removed[sortIdxs, ]

    if (sDiff > 0) {
        removed <- list(removed, NULL)
        rest <- list(normalized, sets[[2]])
    } else if (sDiff < 0) {
        removed <- list(NULL, removed)
        rest <- list(sets[[1]], normalized)
    }  # we already dealt with the improbable equal case

    return (list(removed=removed, rest=rest))
}

.clusterizeShifts <- function (left, right, maxDiff=74)
{
    allShifts <- c(left, right)
    reducedShifts <- reduce(allShifts)
    start(reducedShifts) <- start(reducedShifts) - maxDiff
    end(reducedShifts) <- end(reducedShifts) + maxDiff
    reduce(reducedShifts)
}

.toIRanges <- function (set)
    if (class(set) == "RangedData") {
        IRanges(start=start(set), end=end(set))
    } else {
        set
    }

.dyadPos <- function(ran)
{
    return(round((start(ran) + end(ran)) / 2))
}

.setSizeTo <- function(set, readSize)
{
    return(IRanges(start=.dyadPos(set) - (readSize/2), width=readSize))
}
