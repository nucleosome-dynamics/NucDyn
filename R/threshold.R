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

.fromHss <- function (hss, mc.cores=1)
{
    max.by.rep <- lapply(hss, dlply, "chr", function(x) max(x$end))
    maxs <- do.call(.nmapply, c(max, max.by.rep))
    lines <- lapply(hss, .processHs, maxs, mc.cores=mc.cores)
    lapply(do.call(.averageThreshs, lines), as.data.frame)
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
    VariableThreshold(x=.fromHss(hss, mc.cores=mc.cores))
}

setMethod(
    "applyThreshold",
    signature(threshold="numeric"),
    function (hs, threshold)
        hs[hs$nreads >= threshold, ]
)

setMethod(
    "applyThreshold",
    signature(threshold="integer"),
    function (hs, threshold)
        applyThreshold(hs, as.numeric(threshold))
)

setMethod(
    "applyThreshold",
    signature(threshold="character"),
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
