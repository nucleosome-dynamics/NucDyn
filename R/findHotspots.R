#' Find hotspots in a NucDyn object.
#'
#' Find hotspots from a given nucleosome dynamics.
#'
#' This function is aimed to help in the analysis of the nucleosome dynamics by
#' pointing out these regions with relevant changes in the individual position
#' of the nucleosomes.
#'
#' There are 4 types of basic hotspots, plus their combinations. Basic ones
#' are:
#'
#' * Translational movement of nucleosomes, downstream (+).
#' * Translational movement of nucleosomes, upstream (-).
#' * Nucleosome reads removed from that locus.
#' * Nucleosome reads added to that locus.
#'
#' As translational and coverage changes can happen anywhere, only those
#' involving a certain number of reads are reported. This number can by
#' adjusted by the `threshold` parameter. If `threshold` is a `character`
#' vector representing a percentage value (ie, `"60\\%"`), this will be
#' automatically converted to the absolute value given by the corresponding
#' percentile of the coverage in the window. If, instead, `threshold` is a
#' `numeric` value, this value will be used as absolute threhold.
#'
#' It two adjacent hotspots with shifts in opposite directions are detected but
#' one of them is relatively small in comparison with the other, but will be
#' reported as shifts, disregarding the value of `combined`. We consider two
#' hotspots of the same magnitude if the ratio between the number of reads in
#' one and the other is smaller than `same.magnitude`. This ratio is always
#' performed by using the larger number as numerator and the smaller as
#' denominator; therefore, `same.magnitude` must always be greater of equal
#' than 1.
#' 
#' For example, with `same.magnitude=2`, we consider that 25 reads shifting
#' downstream followed bby 17 reads shifting upstream will be of the same
#' magnitude (25/17 == 1.47 < 2) and we will annotate it as a "DISPERSION". In
#' another example, if we have 25 shifts downstream followed by only 5 shifts
#' upstream (25/5 == 5 > 2), both hotspot will be annotated as "SHIFT".
#'
#' @param dyn NucDyn object with the dynamic to analyze.
#' @param wins Size of the window in base-pairs where the relative scores are
#'   computed
#' @param indel.threshold Maximum p-value for an insertion or delition hotspot
#'   to be considered significant.
#' @param shift.threshold Maximum p-value for shift hotspot to be considered
#'   significant.
#' @param mc.cores If `parallel` support, the number of cores available. This
#'   option is only used if the provided sets are from more than one
#'   chromosome.
#'
#' @return A `data.frame` with the following columns:
#'
#' * chrom: Chromosome name.
#' * coord: Genomic coordinates (average dyad position of affected
#'   nucleosomes).
#' * type: The type of the hotspot (as listed above).
#' * nreads: Number of reads involved in the hotspot.
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
#' @keywords manip
#' @rdname findHotspots
#' @export findHotspots
#'
setGeneric(
    "findHotspots",
    function (dyn, wins=10000, indel.threshold=NULL, shift.threshold=NULL,
              mc.cores=1)
        standardGeneric("findHotspots")
)

#' @rdname findHotspots
#' @importMethodsFrom GenomeInfoDb seqnames
#' @importMethodsFrom IRanges ranges
setMethod(
    "findHotspots",
    signature(dyn="NucDyn"),
    function (dyn, wins=10000, indel.threshold=NULL, shift.threshold=NULL,
              mc.cores=1)
    {
        setA <- set.a(dyn)
        setB <- set.b(dyn)

        chrs <- levels(seqnames(setA)$originals)

        chrIter <- function (chr)
        {
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

       if (length(chrs)) {
           hsLs <- .xlapply(chrs, chrIter)
           hs <- do.call("rbind", hsLs)
       } else {
           hs <- data.frame(start  = integer(),
                            end    = integer(),
                            peak   = integer(),
                            nreads = numeric(),
                            score  = numeric(),
                            type   = character(),
                            chr    = character())
       }

        if (!is.null(indel.threshold) & !is.null(shift.threshold)) {
            hs <- applyThreshold(hs, indel.threshold, shift.threshold)
        }
        hs
    }
)

.calcDiff <- function (x, y)
{
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

.splitBySign <- function (xs)
{
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

.makeVectsEqual <- function(x, y)
{
    # Make two numerical vectors the same size by appending zeros to the
    # shorter one
    vs <- list(x, y)
    lens <- sapply(vs, length)
    d <- max(lens) - min(lens)
    i <- which.min(lens)
    vs[[i]] <- c(vs[[i]], rep(0, d))
    return(vs)
}

#' @importFrom IRanges coverage
.getEqCovs <- function (xs)
    do.call(.makeVectsEqual, lapply(xs, function (x) as.vector(coverage(x))))

.ranScorer <- function (start, end, xs)
    mapply(function (s, e) abs(xs[s:e]), start, end)

.weightedMean <- function (x, weights)
{
    logx <- -log10(x)
    avg <- sum(weights/sum(weights) * logx)
    return (10^(-avg))
}

#' @importMethodsFrom IRanges start end
.ranIter <- function (ran, f, ...)
    mapply(function (i, j) do.call(f, lapply(list(...), `[`, i:j)),
           start(ran),
           end(ran))

.meanArround <- function (x, a, n=3)
{
    i <- (x-3):(x+3)
    i <- i[i > 0]
    vals <- a[i]
    vals <- vals[!is.na(vals)]
    mean(vals)
}

#' @importMethodsFrom IRanges start end
.ran2df <- function (r, xs, pval)
{
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

#' @importFrom plyr rbind.fill
#' @importMethodsFrom nucleR filterFFT
.hsFromCov <- function(x, pvals, names)
{
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

#' @importMethodsFrom IRanges coverage
.getEqualCovs <- function (xs)
    do.call(
        .makeVectsEqual,
        lapply(xs, function (x) as.vector(coverage(x)))
    )

#' Calculate a vector of p-values expressing the difference between two
#' coverages.
#'
#' Calculate a vector of p-values expressing the difference between two
#' coverages. Works by windows.
#'
#' @param x Coverage of the first sample for a given chromosome.
#' @param y Coverage of the second sample for the same chromosome as x.
#' @param wins Size of the window.
#'
#' @return A `numeric` vector of p-values per base-pair.
#'
#' @examples
#'     data(sample_chrII)
#'     pval <- findPVals(sample_chrII[[1]], sample_chrII[[2]], win=10000)
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}, 
#'     Diana Buitrago, Diego Gallego
#' @keywords manip
#' @export findPVals
#'
findPVals <- function (x, y, wins=10000)
    doBySplitting(get_pvals, wins=wins, x, y)

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

#' Apply threshold
#'
#' Apply an indel threshold and a shift threshold to the hotspots
#'
#' @param hs Hotspots returned by findHotspots.
#' @param indel.thresh threshold for the indels
#' @param shift.thresh threshold fo the shifts
#'
#' @return a hotspots `data.frame` filered by the thresholds
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
#' @keywords manip
#'
applyThreshold <- function (hs, indel.thresh, shift.thresh)
{
    is.indel <- hs$type == "EVICTION" | hs$type == "INCLUSION"
    is.shift <- hs$type == "SHIFT +"  | hs$type == "SHIFT -"
    sign.indel <- is.indel & hs$score <= indel.thresh
    sign.shift <- is.shift & hs$score <= indel.thresh
    hs[sign.indel | sign.shift, ]
}
