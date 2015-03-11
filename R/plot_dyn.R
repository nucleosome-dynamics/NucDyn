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
    covA <- IRanges::as.vector(coverage(subA))
    covB <- IRanges::as.vector(coverage(subB)) * norm.factor

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

.parseShifts <- function(shiftList, plot.range)
{
    if (0 %in% lapply(shiftList, length)) {
        return(IRanges())
    } else {
        shiftRange <- start (shiftList[[1]]) >= plot.range[1] &
                      end   (shiftList[[1]]) <= plot.range[2]
        shiftList <- lapply(shiftList, function(set) set[shiftRange])
        shiftPos <- .dyadPos(shiftList[[1]])
        shiftDiff <- .dyadPos(shiftList[[2]]) - shiftPos

        shiftRan <- IRanges(start=shiftPos, end=(shiftPos + abs(shiftDiff)))

        return(shiftRan)
    }
}

.subsetCov <- function(covRan, plot.range)
{
    if (is.null(covRan)) {
        return(IRanges())
    } else {
        subset <- start (covRan) >= plot.range[1] &
                  end   (covRan) <= plot.range[2]
        covRan <- covRan[subset]

        return(covRan)
    }
}

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
        covDel <- IRanges::as.vector(coverage(deletions))
        lines(covDel, type="h", lwd=2, col=args[["col.indel"]][2])
    }

    if (length(insertions)) {
        covIns <- IRanges::as.vector(coverage(insertions))
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
         c("", "", "Inserted reads", "Removed reads", "Upstream shift", "Downstream shift"),
         col=c(NA, NA, args[["col.indel"]], args[["col.arr"]]),
         lty=c(NA, NA, 1, 1, 1, 1), lwd=c(NA, NA, 3, 3, 3, 3), bty="n", cex=0.8)
}


plotDynamics <- function(dyn, plot.range=NULL, chr=NULL,
                         dyn.name="Dyn", expA.name="Ref 1", expB.name="Ref 2",
                         norm.factor=1, ...)
invisible({
    if (is.null(chr)) {
        warning("No chromosome specified. If the input contained more than one chromosome, they will apprear concatenated in the plot.")
        dynRanges <- lapply(list(set.a(dyn), set.b(dyn)), ranges)
    } else {
        dynRanges <- lapply(
            list(set.a(dyn), set.b(dyn)),
            function (x) ranges(x[seqnames(x) == chr])
        )
    }

    types <- unique(c(names(dynRanges[[1]]), names(dynRanges[[2]])))
    dyn <- lapply(
        types,
        function (t)
            lapply(dynRanges, function(x) x[[t]])
    )
    names(dyn) <- types

    if (is.null(plot.range)) {
        ran <- range(do.call(c, lapply(dynRanges, GenomicRanges::unlist)))
        plot.range <- c(start(ran), end(ran))
    }

    .plotCoverage(dyn$originals, plot.range=plot.range,
                  dyn.name=dyn.name, expA.name=expA.name, expB.name=expB.name,
                  norm.factor=norm.factor, ...)

    .addDynamics(dyn, plot.range=plot.range, ...)
})
