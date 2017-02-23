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
    round((start(ran) + end(ran)) / 2)

.setSizeTo <- function(set, readSize)
    IRanges(start=.dyadPos(set) - (readSize/2), width=readSize)

.initEmptyVect <- function(n)
    as.integer(rep(0, n))

.dlplyf <- function (data, splitter, fun, ...)
    # Similar to dlply, but splits the data using a splitter function that
    # should return a list of data.frames
    lapply(splitter(data), fun, ...)

.ddplyf <- function (data, splitter, fun, ...)
    # Similar to ddply, but splits the data using a splitter function that
    # should return a list of data.frames
    do.call(rbind, .dlplyf(data, splitter, fun, ...))

.xdlply_rep <- function (X, VAR, FUN, ..., report=FALSE, mc.cores=1)
{   # multicore version of dlply that also reports the name of the element
    # being processed
    xs <- unique(X[[VAR]])
    res <- .xlapply(xs,
                    function (i) {
                        if (report) {
                            message("Starting ", i)
                        }
                        res <- FUN(X[X[[VAR]] == i, ], ...)
                        if (report) {
                            message(i, " done")
                        }
                        res
                    },
                    mc.cores=mc.cores)
    names(res) <- xs
    res
}

.xddply_rep <- function (X, VAR, FUN, ..., report=FALSE, mc.cores=1)
    # multicore version of ddply that also reports the name of the element
    # being processed
    do.call(rbind,
            .xdlply_rep(X,
                        VAR,
                        FUN,
                        report=report,
                        ...,
                        mc.cores=mc.cores))

.nmapply <- function (FUN, ..., MoreArgs=NULL)
{   # Similar to mapply but using the names of the elements
    args <- list(...)
    ns <- unique(unlist(lapply(args, names)))
    res <- lapply(ns,
                  function (n)
                      do.call(FUN,
                              c(lapply(args, `[[`, n),
                                MoreArgs)))
    names(res) <- ns
    res
}

.vectorMean <- function (...)
{   # Get the mean of a number of numeric vectors
    args <- list(...)
    Reduce(`+`, args) / length(args)
}

compose <- function(...)
{   # Function composition.
    comp2 <- function(f, g) {
        force(f)
        force(g)
        function(...) f(g(...))
    }
    Reduce(comp2, list(...))
}
