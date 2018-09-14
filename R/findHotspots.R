#' Find hotspots in a NucDyn object.
#'
#' Find hotspots from a given nucleosome dynamics.
#'
#' This function is aimed to help in the analysis of the nucleosome dynamics by
#' pointing out these regions with relevant changes in the individual position
#' of the nucleosomes.
#'
#' There are 4 types of basic hotspots:
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
#' @param nuc list of two `GRanges` objects containing the nucleosome calls of 
#'   experiment 1 and experiment 2. 
#' @param wins Size of the window in base-pairs where the relative scores are
#'   computed
#' @param indel.threshold Maximum p-value for an `INCLUSION` or `EVICTION` 
#'   hotspot to be considered significant.
#' @param shift.threshold Maximum p-value for a `SHIFT +` or `SHIFT -`
#'   hotspot to be considered significant.
#' @param indel.nreads Minimum number of reads in an `INCLUSION` or `EVICTION` 
#'   hotspot.
#' @param shift.nreads Minimum number of reads in a `SHIFT +` or `SHIFT -`
#'   hotspot.
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
#' @examples
#'     data(readsG2_chrII)
#'     data(readsM_chrII)
#'     data(nuc_chrII)
#'     dyn <- nucleosomeDynamics(setA=readsG2_chrII, setB=readsM_chrII)
#'     findHotspots(dyn, nuc_chrII)
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}, 
#'     Diana Buitrago \email{diana.buitrago@@irbbarcelona.org}
#' @keywords manip
#' @rdname findHotspots
#' @export findHotspots
#'
setGeneric(
    "findHotspots",
    function (dyn, nuc, wins=10000, indel.threshold=NULL, shift.threshold=NULL,
              indel.nreads=NULL, shift.nreads=NULL, mc.cores=1)
        standardGeneric("findHotspots")
)

#' @rdname findHotspots
#' @importMethodsFrom GenomeInfoDb seqnames
#' @importMethodsFrom IRanges ranges
#' @importFrom GenomicRanges findOverlaps
setMethod(
    "findHotspots",
    signature(dyn="NucDyn"),
    function (dyn, nuc, wins=10000, indel.threshold=NULL, shift.threshold=NULL,
              indel.nreads=NULL, shift.nreads=NULL, mc.cores=1)
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

        hs_gr <- GRanges(hs$chr,IRanges(start=hs$start, end=hs$end), 
                         peak=hs$peak, nreads=hs$nreads, score=hs$score, 
                         type=hs$type)

        hs_shiftp <- hs_gr[grep("SHIFT \\+", hs_gr$type)] 
        hs_shiftm <- hs_gr[grep("SHIFT -", hs_gr$type)] 

        #Project shift to dyad
        nuc1 <- nuc[[1]]

        ovlp_nuc1p <- as.data.frame(findOverlaps(nuc1, hs_shiftp))       
        ovlp_nuc1p$start <- start(nuc1[ovlp_nuc1p$queryHits])
        ovlp_nuc1p <- lapply(split(ovlp_nuc1p, ovlp_nuc1p$subjectHits),
           function(x) {
              tmp_nuc = nuc1[x$queryHits[which.min(x$start)]]
              if(hs_shiftp[x$subjectHits[1]]$type=="SHIFT +1"){
                 tmp_nuc = nuc1[x$queryHits[which.max(x$start)]]
                 tmp_shf = hs_shiftp[x$subjectHits[which.max(x$start)]]
              } else {
                 tmp_nuc = nuc1[x$queryHits[which.min(x$start)]]
                 tmp_shf = hs_shiftp[x$subjectHits[which.min(x$start)]]
              }
              data.frame(chr    = as.character(seqnames(tmp_nuc)),
                         start  = as.numeric(start(tmp_nuc)),
                         end    = as.numeric(end(tmp_nuc)),
                         peak   = tmp_shf$peak, 
                         nreads = tmp_shf$nreads,
                         score  = tmp_shf$score,
                         type   = "SHIFT +") #tmp_shf$type)
           })

        hs_shiftp_nuc <- do.call(rbind, ovlp_nuc1p)
                   
        ovlp_nuc1m <- as.data.frame(findOverlaps(nuc1, hs_shiftm))        
        ovlp_nuc1m$end <- end(nuc1[ovlp_nuc1m$queryHits])
        ovlp_nuc1m <- lapply(split(ovlp_nuc1m, ovlp_nuc1m$subjectHits),
           function(x) {
              if(hs_shiftm[x$subjectHits[1]]$type=="SHIFT -1"){
                 tmp_nuc = nuc1[x$queryHits[which.min(x$end)]]
                 tmp_shf = hs_shiftm[x$subjectHits[which.min(x$end)]]
              } else {
                 tmp_nuc = nuc1[x$queryHits[which.max(x$end)]]
                 tmp_shf = hs_shiftm[x$subjectHits[which.max(x$end)]]
              }
              data.frame(chr    = as.character(seqnames(tmp_nuc)),
                         start  = as.numeric(start(tmp_nuc)),
                         end    = as.numeric(end(tmp_nuc)),
                         peak   = tmp_shf$peak,
                         nreads = tmp_shf$nreads,
                         score  = tmp_shf$score,
                         type   = "SHIFT -") #tmp_shf$type)
           })
        
        hs_shiftm_nuc <- do.call(rbind, ovlp_nuc1m)
        
        #append indels
        hs_in <- hs[hs$type=="INCLUSION",]
        hs_del <- hs[hs$type=="EVICTION",]

        hs_out <- rbind(hs_in, hs_del, hs_shiftp_nuc, hs_shiftm_nuc)                

        #apply thresholds
        if (!is.null(indel.threshold) & !is.null(shift.threshold) &
            !is.null(indel.nreads) & !is.null(shift.nreads) ) {
            hs_out <- applyThreshold(hs_out, indel.thresh=indel.threshold, 
                                     shift.thresh=shift.threshold,
                                     indel.nreads=indel.nreads, 
                                     shift.nreads=shift.nreads)
        }
        #remove overlappyng shift + and -
        id_hs <- paste(hs_out$chr, hs_out$start, hs_out$end, sep="_")
        is.sp <- hs_out$type == "SHIFT +"
        is.sm <- hs_out$type == "SHIFT -"
        
        rm_sp <- id_hs[is.sp][id_hs[is.sp]%in%id_hs[is.sm]]
        rm_sm <- id_hs[is.sm][id_hs[is.sm]%in%id_hs[is.sp]]
        rm_shift <- c(which(id_hs %in% rm_sp), which(id_hs %in% rm_sm))
        if(length(rm_shift)>0){
           hs_out <- hs_out[-c(which(id_hs%in%rm_sp), which(id_hs%in%rm_sm)),]
        }
        return(hs_out)
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
#'     Diana Buitrago \email{diana.buitrago@@irbbarcelona.org}, 
#'     Diego Gallego
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
                               dyn[["right.shifts"]][[2]]))
    diff <- do.call(`-`, covs)
    shiftsp <- .hsFromCov(diff, pvals, c("SHIFT +1", "SHIFT +2"))

    ###########################################################################

    covs <- .getEqualCovs(list(dyn[["left.shifts"]][[1]],
                               dyn[["left.shifts"]][[2]]))
    diff <- do.call(`-`, covs)
    shiftsm <- .hsFromCov(diff, pvals, c("SHIFT -1", "SHIFT -2"))

    ###########################################################################

    hs <- rbind(indels, shiftsp, shiftsm)
    hs[order(hs$start), ]
}

#' Apply threshold
#'
#' Apply an indel threshold and a shift threshold to the hotspots
#'
#' @param hs Hotspots returned by findHotspots.
#' @param indel.thresh threshold for the indels.
#' @param shift.thresh threshold for the shifts.
#' @param indel.nreads minimum number of reads for the indels.
#' @param shift.nreads minimum number of reads for the shifts.
#'
#' @return a hotspots `data.frame` filered by the thresholds.
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}, 
#'     Diana Buitrago \email{diana.buitrago@@irbbarcelona.org}
#' @keywords manip
#'
applyThreshold <- function(hs, indel.thresh, shift.thresh, 
                           indel.nreads, shift.nreads)
{
    is.indel <- hs$type == "EVICTION" | hs$type == "INCLUSION"
    is.shift <- hs$type == "SHIFT +"  | hs$type == "SHIFT -"
    sign.indel <- is.indel & hs$score <= indel.thresh & hs$nreads >= indel.nreads
    sign.shift <- is.shift & hs$score <= shift.thresh & hs$nreads >= shift.nreads
    hs[sign.indel | sign.shift, ]
}
