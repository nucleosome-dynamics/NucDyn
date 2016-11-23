.selfZip <- function (x, f, ...)
    mapply(f,
           x[-length(x)],
           x[-1],
           MoreArgs=list(...),
           SIMPLIFY=FALSE)

.islandSplitDiscrete <- function (x) {
    if (sum(x) == 0) {
        return(c())
    } else {
        nums <- c(x, 0, 0)
        sign.diffs <- diff(sign(diff(nums)))

        a <- sign.diffs[-length(sign.diffs)]
        b <- sign.diffs[-1]

        all.vals <- ifelse(a > 0 & b < 0,
                           "up",
                           ifelse(a < 0 & b > 0,
                                  "down",
                                  "_"))
        sel <- all.vals != "_"
        pos <- which(sel)
        val <- all.vals[sel]

        val.a <- val[-length(val)]
        val.b <- val[-1]
        pos.a <- pos[-length(pos)]
        pos.b <- pos[-1]

        middles <- ifelse(val.a == "down" & val.b == "up",
                          round((pos.a + pos.b)/2),
                          0)

        middles[middles != 0]
    }
}

.islandSplitContinuous <- function (x) {
    if (sum(x) == 0) {
        return(c())
    } else {
        nums <- c(x, 0, 0)
        sign.diffs <- diff(sign(diff(nums)))
        return(which(sign.diffs == 2))
    }
}

.getIslandLims <- function (xs) {
    nums <- c(xs, 0, 0)
    sign.diffs <- diff(sign(diff(nums)))
    is.zero <- nums[2:(length(nums) - 1)] == 0
    sel <- sign.diffs > 0 & is.zero
    return(which(sel))
}

.getIslandedLims <- function (xs) {
    doRange <- function (i, j, xs) {
        island <- xs[(i+2):j]
        disc.splits <- .islandSplitDiscrete(island)
        cont.splits <- .islandSplitContinuous(island)
        c(disc.splits, cont.splits) + i
    }

    extraIsland.lims <- .getIslandLims(xs)
    intraIsland.lims <- unlist(.selfZip(extraIsland.lims, doRange, xs))

    sort(unique(c(extraIsland.lims, intraIsland.lims)))
}

.lims2range <- function (lims)
    IRanges(start=c(1, lims[-length(lims)] + 2), end=lims)

.peaksFromRanges <- function (xs, rans) {
    doRange <- function (start, end, xs) {
        sub.ran <- xs[start:end]
        if (sum(sub.ran) == 0) {
            return(NULL)
        } else {
            max.vals <- which(sub.ran == max(sub.ran))
            center.pos <- round(mean(max.vals))
            return(center.pos + start)
        }
    }
    unlist(mapply(doRange,
                  start(rans),
                  end(rans),
                  MoreArgs=list(xs),
                  SIMPLIFY=FALSE))
}

.filterEmptyRanges <- function (rans, xs) {
    doRange <- function (start, end, xs)
        sum(xs[start:end]) > 0
    sel <- mapply(doRange,
                  start(rans),
                  end(rans),
                  MoreArgs=list(xs))
    rans[sel]
}

.getHsRanges <- function (x)
    .filterEmptyRanges(.lims2range(.getIslandedLims(x)), x)

