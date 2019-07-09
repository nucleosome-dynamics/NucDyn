#' Plot a NucDyn object.
#'
#' Plot a visual representation of a NucDyn object.
#' Coverage profile for ref1 will be shown as a solid grey background and
#' coverage profile for ref2 will be shown as a dotted profile.  Superimposed,
#' arrows showing read shifts (upstream in blue and downstream in red). Indels
#' will appear as a small coverage profile at the bottom part of the plot
#' (insertions in green and deletions in red).
#'
#' @param dyn `NucDyn` object with the dynamic to plot.
#' @param plot.range Range from the `NucDyn` object to plot. If not specified,
#'     the whole set will be plotted. If `dyn` contains more than one
#'     chromosome, they will appear concatenated in the plot.
#' @param chr Chromosome of the `NucDyn` object to plot. If not specified,
#'     all chromosomes will appear plotted concatenated.
#' @param dyn.name Name to be given to the dyanamics that will be displayed in
#'     the plot.
#' @param expA.name Name to be given to the first data set of the dyanamics
#'     that will be displayed in the plot.
#' @param expB.name Name to be given to the second data set of the dyanamics
#'     that will be displayed in the plot.
#' @param norm.factor Normalization factor between ref1 and ref2. Use it to
#'     visualize both coverages profiles on a similar scale if one of them has
#'     a significantly higher coverage.
#' @param \dots Other parameters passed to [graphics::plot()] function.
#'
#' @return Void
#'
#' @author Oscar Flores, Ricard Illa
#'     Diana Buitrago \email{diana.buitrago@@irbbarcelona.org}
#' @rdname plotDynamics
#' @export plotDynamics
#' @keywords hplot
#
setGeneric(
    "plotDynamics",
    function(dyn, ...) standardGeneric("plotDynamics")
)

#' @rdname plotDynamics
#' @importMethodsFrom IRanges ranges start end
#' @importMethodsFrom GenomeInfoDb seqnames
setMethod(
    "plotDynamics",
    signature(dyn="NucDyn"),
    function(dyn, plot.range=NULL, chr=NULL,
             dyn.name="Dyn", expA.name="Ref 1", expB.name="Ref 2",
             norm.factor=1, ...) {
        invisible({
            if (is.null(chr)) {
                warning("No chromosome specified. ", 
                        "If the input contained more than one chromosome, ",
                        "they will apprear concatenated in the plot.")
                dynRanges <- lapply(list(set.a(dyn), set.b(dyn)), ranges)
            } else {
                dynRanges <- lapply(
                    list(set.a(dyn), set.b(dyn)),
                    function (x) ranges(x[seqnames(x) == chr])
                )
            }

            types <- unique(c(names(dynRanges[[1]]),
                              names(dynRanges[[2]])))
            dyn <- lapply(
                types,
                function (t)
                    lapply(dynRanges, function(x) x[[t]])
            )
            names(dyn) <- types

            if (is.null(plot.range)) {
                ran <- range(do.call(c,
                                     lapply(dynRanges,
                                            unlist)))
                plot.range <- c(start(ran), end(ran))
            }

            .plotCoverage(dyn$originals, plot.range=plot.range,
                          dyn.name=dyn.name, expA.name=expA.name,
                          expB.name=expB.name, norm.factor=norm.factor, ...)

            .addDynamics(dyn, plot.range=plot.range, ...)
        })
    }
)

#' @importMethodsFrom IRanges start end coverage
#' @importFrom graphics plot lines legend
.plotCoverage <- function(originalReads, plot.range,
                          dyn.name="Dyn", expA.name="Ref 1", expB.name="Ref 2",
                          norm.factor=1, ...)
{
    args <- list(...)

    expA <- originalReads[[1]]
    expB <- originalReads[[2]]

    name <- paste0(dyn.name, " (", expA.name, " -> ", expB.name, ")")

    subA <- expA[start(expA) >= plot.range[1] & end(expA) <= plot.range[2]]
    subB <- expB[start(expB) >= plot.range[1] & end(expB) <= plot.range[2]]

    if (!length(subA) || !length(subB)) {
        stop("No reads present in the chromosome/range asked")
    }

    # set defaults for unset parameters
    if (is.null(args[["xlim"]])) {
        args[["xlim"]] <- c(plot.range[1], plot.range[2])
    }
    if (is.null(args[["main"]])) {
        args[["main"]] <- name
    }
    if (is.null(args[["xlab"]])) {
        args[["xlab"]] <- "position"
    }
    if (is.null(args[["ylab"]])) {
        args[["ylab"]] <- "coverage"
    }
    if (is.null(args[["col"]])) {
        args[["col"]] <- c("grey", "black")
    } else if (length(args[["col"]]) < 2) {
        args[["col"]] <- rep(args[["col"]], 2)
    }

    # calculate coverages
    covA <- as.vector(coverage(subA))
    covB <- as.vector(coverage(subB)) * norm.factor

    args1 <- args
    # remove composite pars
    args1 <- args1[grep(".", names(args1), fixed=TRUE, invert=TRUE)]
    args1[["col"]] <- args[["col"]][1]
    args1[["x"]] <- covA
    args1[["type"]] <- "h"
    args1[["lty"]] <- 1
    args1[["lwd"]] <- 2

    do.call(plot, args1)
    lines(covB, lwd=2, col=args[["col"]][2], lty=3)

    legend("topleft",
           c(paste("Ref1:", expA.name), paste("Ref2:", expB.name)),
           col=args[["col"]][1:2], lty=c(1,3), lwd=c(3,3), bty="n", cex=0.8)
}

#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges start end
.parseShifts <- function(shiftList, plot.range)
{
    if (0 %in% lapply(shiftList, length)) {
        return(IRanges())
    } else {
        shiftRange <- start(shiftList[[1]]) >= plot.range[1] &
                      end  (shiftList[[1]]) <= plot.range[2]
        shiftList <- lapply(shiftList, function(set) set[shiftRange])
        shiftPos <- .dyadPos(shiftList[[1]])
        shiftDiff <- .dyadPos(shiftList[[2]]) - shiftPos

        shiftRan <- IRanges(start=shiftPos, end=(shiftPos + abs(shiftDiff)))

        return(shiftRan)
    }
}

#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges start end
.subsetCov <- function(covRan, plot.range)
{
    if (is.null(covRan)) {
        return(IRanges())
    } else {
        subset <- start (covRan) >= plot.range[1] &
                  end(covRan) <= plot.range[2]
        covRan <- covRan[subset]

        return(covRan)
    }
}

#' @importFrom graphics par lines arrows
#' @importMethodsFrom IRanges coverage disjointBins start end
.addDynamics <- function(dyn, plot.range, ...)
{
    args <- list(...)

    left.shifts <- .parseShifts(dyn$left.shifts, plot.range)
    right.shifts <- .parseShifts(dyn$right.shifts, plot.range)

    # set defaults for unset parameters
    if (is.null(args[["col"]])) {
        args[["col"]] <- c("grey", "black")
    } else if (length(args[["col"]]) < 2) {
        args[["col"]] <- rep(args[["col"]], 2)
    }
    if (is.null(args[["col.arr"]])) {  # arrow colors
        args[["col.arr"]] <- c("darkred", "darkblue")
    } else if (length(args[["col.arr"]]) < 2) {
        args[["col.arr"]] <- rep(args[["col.arr"]], 2)
    }
    if (is.null(args[["col.indel"]])) {  # indel colors
        args[["col.indel"]] <- c("green", "red")
    } else if (length(args[["col.indel"]]) < 2) {
        args[["col.indel"]] <- rep(args[["col.indel"]], 2)
    }

    if (is.null(args[["arr.length"]])) {
        args[["arr.length"]] <- 0.05
    }
    if (is.null(args[["arr.lwd"]])) {
        args[["arr.lwd"]] <- 1
    }
    if (is.null(args[["arr.angle"]])) {
        args[["arr.angle"]] <- 30
    }
    if (is.null(args[["arr.code"]])) {
        args[["arr.code"]] <- 2
    }

    par.usr <- par("usr")
    par.ypc <- (par.usr[4] - par.usr[3]) / 100

    deletions <- .subsetCov(dyn$indels[[1]], plot.range)
    insertions <- .subsetCov(dyn$indels[[2]], plot.range)

    if (length(deletions) > 0) {
        covDel <- as.vector(coverage(deletions))
        lines(covDel, type="h", lwd=2, col=args[["col.indel"]][2])
    }

    if (length(insertions)) {
        covIns <- as.vector(coverage(insertions))
        lines(covIns, type="h", lwd=2, col=args[["col.indel"]][1])
    }

    middle <- (par.usr[4] - par.usr[3]) * 0.5

    db.left <- disjointBins(left.shifts)
    db.right <- disjointBins(right.shifts)

    for (i in 1:length(right.shifts)) {
        x1 <- start(right.shifts)[i]
        y1 <- middle + par.ypc*db.right[i]
        x2 <- end(right.shifts)[i]
        y2 <- middle + par.ypc*db.right[i]

        arrows(x1, y1, x2, y2,
               length=args[["arr.length"]], col=args[["col.arr"]][2],
               angle=args[["arr.angle"]],   code=args[["arr.code"]],
               lwd=args[["arr.lwd"]])
    }

    for (i in 1:length(left.shifts)) {
        x1 <- end(left.shifts)[i]
        y1 <- middle - par.ypc*db.left[i]
        x2 <- start(left.shifts)[i]
        y2 <- middle - par.ypc*db.left[i]

        arrows(x1, y1, x2, y2,
               length=args[["arr.length"]], col=args[["col.arr"]][1],
               angle=args[["arr.angle"]],   code=args[["arr.code"]],
               lwd=args[["arr.lwd"]])
    }

    legend("topleft",
         c("", "", "Inserted reads", "Removed reads", "Upstream shift",
           "Downstream shift"),
         col=c(NA, NA, args[["col.indel"]], args[["col.arr"]]),
         lty=c(NA, NA, 1, 1, 1, 1), lwd=c(NA, NA, 3, 3, 3, 3), bty="n",
         cex=0.8)
}
