.dfToLine <- function(df, max.vals, chr=NULL, expand=30) {
    cWrap <- function (vals, max.val, expand=0) {
        start <- vals$start - expand
        start[start < 1] <- 1
        start <- as.integer(start)

        end <- vals$end + expand
        end[end > max.val] <- max.val
        end <- as.integer(end)

        xs <- as.integer(vals$nreads)
        n <- as.integer(length(xs))
        out <- as.numeric(rep(0, max.val))

        cOut <- .C("cov_ranges",
                   xs, start, end, n,
                   max.val, out=out)$out
        return(cOut)
    }

    if (is.null(chr)) {
        chr <- df[1, "chr"]
    }
    max.val <- max.vals[[chr]]

    nonzero <- df[df$nreads != 0, ]
    xs <- cWrap(nonzero, max.val, expand)
    filterFFT(xs,
              useOptim=TRUE,
              pcKeepComp=0.005)
}

.averageThreshs <- function (...)
    .nmapply(function (...) .nmapply(.vectorMean, ...), ...)

.processHs <- function (hs, maxs, mc.cores=1)
    .xdlply_rep(hs,
                "chr",
                function (x)
                    .dlplyf(x,
                            .typeSplitter,
                            .dfToLine,
                            maxs),
                mc.cores=mc.cores)

variableThreshFromHss <- function (..., mc.cores=1)
{
    hss <- list(...)
    max.by.rep <- lapply(hss, dlply, "chr", function(x) max(x$end))
    maxs <- do.call(.nmapply, c(max, max.by.rep))
    lines <- lapply(hss, .processHs, maxs, mc.cores=mc.cores)
    t <- lapply(do.call(.averageThreshs, lines), as.data.frame)
    VariableThreshold(x=t, scale=3)
}

getVariableThreshold <- function (..., mc.cores=1)
{   # Expects pairs of replicates
    message("Calculating dynamics")
    dyns <- lapply(list(...),
                   function (x)
                       nucleosomeDynamics(x[[1]],
                                          x[[2]],
                                          mc.cores=mc.cores))
    message("Finding hotspots")
    hss <- lapply(dyns,
                  findHotspots,
                  combined=FALSE,
                  threshold=NULL,
                  mc.cores=mc.cores)
    message("Averaging the different replicates")
    do.call(variableThreshFromHss,
            c(hss, list(mc.cores=mc.cores)))
}

setMethod(
    "applyThreshold",
    signature(hs="data.frame", threshold="numeric"),
    function (hs, threshold, scale=NULL)
        hs[hs$nreads >= threshold, ]
)

setMethod(
    "applyThreshold",
    signature(hs="data.frame", threshold="integer"),
    function (hs, threshold)
        applyThreshold(hs, as.numeric(threshold))
)

setMethod(
    "applyThreshold",
    signature(hs="data.frame", threshold="character"),
    function (hs, threshold) {
        relative <- as.numeric(sub("%", "", threshold)) / 100
        ddply(hs,
              "chr",
              .ddplyf,
              .typeSplitter,
              function(df) {
                  thresh <- quantile(df$nreads, relative)
                  df[df$nreads >= thresh, ]
              })
    }
)

setMethod(
    "applyThreshold",
    signature(hs="data.frame", threshold="VariableThreshold"),
    function (hs, threshold, scale=NULL) {

        if (is.null(scale)) {
            scale <- threshold@scale
        }

        applyThreshVect <- function (hs, t) {
            diff <- max(hs$coord) - length(t)
            if (diff > 0) {
                t <- c(t, rep(0, diff))
            }
            hs[hs$nreads >= (t[hs$coord] * scale), ]
        }

        applyAtChr <- function (chr.hs, chr.thresh) {
            res <- do.call(rbind,
                           .nmapply(applyThreshVect,
                                    .typeSplitter(chr.hs),
                                    chr.thresh))
            res[order(res$coord), ]
        }

        res <- do.call(rbind,
                       .nmapply(applyAtChr,
                                dlply(hs, "chr", identity),
                                threshold@x))
        rownames(res) <- NULL
        res
    }
)
