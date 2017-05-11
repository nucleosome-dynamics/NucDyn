.getPVals <- function (x, y) {
    xs <- as.integer(x)
    ys <- as.integer(y)
    n <- as.integer(length(xs))
    out <- as.numeric(rep(0, n))
    cOut <- .C("get_pvals", xs, ys, n, out=out)
    return(cOut$out)
}

.calcDiff <- function (x, y) {
    X <- sum(x)
    Y <- sum(y)

    N <- X + Y
    n <- x + y

    E <- n * (X/N)
    V <- E * (Y/N) * ((N-n)/(N-1))

    z <- (x-E) / sqrt(V)
    z[is.nan(z)] <- 0
    z
}

.splitBySign <- function (xs) {
    # Convert a numeric vector into a list of two numeric vectors.
    # The first one will have the originaly positive numbers set to zero and
    # the second one, the originaly negative ones.
    a.prof <- xs
    a.prof[a.prof < 0] <- 0

    b.prof <- xs
    b.prof[b.prof > 0] <- 0
    b.prof <- abs(b.prof)

    list(a=a.prof, b=b.prof)
}

.catcher <- function (f)
    # Decorator for functions that may fail
    function (...)
        tryCatch(f(...),
                 error=function (e) NULL)

.makeVectsEqual <- function(x, y) {
    # Make two numerical vectors the same size by appending zeros to the
    # shorter one
    vs <- list(x, y)
    lens <- sapply(vs, length)
    d <- max(lens) - min(lens)
    i <- which.min(lens)
    vs[[i]] <- c(vs[[i]], rep(0, d))
    return(vs)
}

.getEqCovs <- function (xs)
    do.call(.makeVectsEqual, lapply(xs, function (x) as.vector(coverage(x))))

.ranScorer <- function (start, end, xs)
    mapply(function (s, e) abs(xs[s:e]), start, end)

.weightedMean <- function (x, weights) {
    logx <- -log10(x)
    avg <- sum(weights/sum(weights) * logx)
    return (10^(-avg))
}

.ranIter <- function (ran, f, ...)
    mapply(function (i, j)
               do.call(f, lapply(list(...), `[`, i:j)),
           start(ran),
           end(ran))

.meanArround <- function (x, a, n=3) {
    i <- (x-3):(x+3)
    i <- i[i > 0]
    vals <- a[i]
    vals <- vals[!is.na(vals)]
    mean(vals)
}

.ran2df <- function (r, xs, pval) {
    if (length(r)) {
        peak <- start(r) + .ranIter(r, which.max, xs)
        score <- .ranIter(r, .weightedMean, pval, xs)
        nreads <- xs[peak]

        score[is.na(score)] <- 1
        nreads[is.na(nreads)] <- 0

        data.frame(start  = start(r),
                   end    = end(r),
                   peak   = peak,
                   nreads = nreads,
                   score  = score)
    } else {
        data.frame(start  = integer(),
                   end    = integer(),
                   peak   = numeric(),
                   nreads = numeric(),
                   score  = numeric())
    }
}

.hsFromCov <- function(x, pvals, names) {
    by.sign <- .splitBySign(x)
    filtered <- lapply(by.sign,
                       .catcher(filterFFT),
                       pcKeepComp=0.01,
                       useOptim=TRUE)
    rans <- lapply(filtered, .getHsRanges)
    dfs <- mapply(.ran2df,
                  rans,
                  filtered,
                  MoreArgs=list(pvals),
                  SIMPLIFY=FALSE)
    for (i in seq_along(dfs)) {
        if (nrow(dfs[[i]])) {
            dfs[[i]][["type"]] <- names[[i]]
        } else {
            dfs[[i]][["type"]] <- character()
        }
    }
    rbind.fill(dfs)
}

.getEqualCovs <- function (xs)
    do.call(.makeVectsEqual,
            lapply(xs,
                   function (x)
                       as.vector(coverage(x))))

findPVals <- function (x, y, wins=10000)
    doBySplitting(.getPVals, wins=wins, x, y)

.findInRange <- function (dyn, wins=10000)
{
    full.covs <- .getEqualCovs(dyn$originals)
    pvals <- findPVals(full.covs[[1]], full.covs[[2]], wins)

    ###########################################################################

    z <- doBySplitting(.calcDiff,
                       wins=wins,
                       full.covs[[1]],
                       full.covs[[2]])

    indels <- .hsFromCov(z, pvals, c("EVICTION", "INCLUSION"))

    ###########################################################################

    covs <- .getEqualCovs(list(dyn[["right.shifts"]][[1]],
                               dyn[["left.shifts"]][[1]]))
    diff <- do.call(`-`, covs)
    shifts <- .hsFromCov(diff, pvals, c("SHIFT +", "SHIFT -"))

    ###########################################################################

    hs <- rbind(indels, shifts)
    hs[order(hs$start), ]
}

applyThreshold <- function (hs, indel.thresh, shift.thresh) {
    is.indel <- hs$type == "EVICTION" | hs$type == "INCLUSION"
    is.shift <- hs$type == "SHIFT +"  | hs$type == "SHIFT -"
    sign.indel <- is.indel & hs$score <= indel.thresh
    sign.shift <- is.shift & hs$score <= indel.thresh
    hs[sign.indel | sign.shift, ]
}

setMethod(
    "findHotspots",
    signature(dyn="NucDyn"),
    function (dyn, wins=10000, indel.threshold=NULL, shift.threshold=NULL, mc.cores=1) {
        setA <- set.a(dyn)
        setB <- set.b(dyn)

        chrs <- levels(seqnames(setA)$originals)

        chrIter <- function(chr, range=NULL) {
            message(paste("Starting", chr))
            f <- function(x) ranges(x[seqnames(x) == chr])
            chrDyn <- mapply(list, f(setA), f(setB), SIMPLIFY=FALSE)
            hs <- .findInRange(chrDyn, wins=wins)
            if (nrow(hs)) {
                hs$chr <- chr
            }
            rownames(hs) <- NULL
            message(paste(chr, "done"))
            hs
        }

        hsLs <- .xlapply(chrs, chrIter, mc.cores=mc.cores)
        hs <- do.call("rbind", hsLs)
        if (!is.null(indel.threshold) & !is.null(shift.threshold)) {
            hs <- applyThreshold(hs, indel.threshold, shift.threshold)
        }
        hs
    }
)
