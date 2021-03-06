#' @importFrom IRanges IRanges
.buildZoneRange <- function(start, end, width)
{
    if (any(sapply(list(start, end), length) == 0)) {
        return(NULL)
    }
    npos <- ceiling((end - start) / width)
    zone.starts <- 0:(npos-1) * width + start
    zone.rans <- IRanges(start=zone.starts, width=width)
    return(zone.rans)
}

#' @importMethodsFrom IRanges start end width "start<-" "end<-"
.shrinker <- function(rans, by=2, except.start=NULL, except.end=NULL)
{
    old.width <- width(rans)[1]
    new.width <- old.width/2

    x <- old.width/2 - new.width/2
    a <- ceiling(x)
    b <- floor(x)

    if (is.null(except.start)) {
        start(rans) <- start(rans) + a
    } else {
        start(rans[-except.start]) <- start(rans[-except.start]) + a
    }

    if (is.null(except.end)) {
        end(rans) <- end(rans) - b
    } else {
        end(rans[-except.end]) <- end(rans[-except.end]) - b
    }

    return(rans)
}

.absPos <- function(sh, sub.rans, sets)
{
    xs <- sh[[1]][as.logical(sh[[1]])]
    f <- function(a, b, c) find_abs_pos(xs, a, b, c)
    mapply(f, sh, sub.rans, sets, SIMPLIFY=FALSE)
}

#' @importMethodsFrom IRanges start end
.getInRange <- function(dyads, sets, ran)
{
    f <- function(a, xs)
        end(a) > start(xs) & start(a) < end(xs)

    pos <- do.call(`&`, lapply(mapply(`[`, sets, dyads), f, ran))
    lapply(dyads, `[`, pos)
}

#' @importMethodsFrom IRanges start end
.doZone <- function(i, bigzones, smallzones, sets, max.dist, min.dist)
{
    f <- function(a, xs)
        end(a) > start(xs) & start(a) < end(xs)

    sub.rans <- lapply(sets, function(s) s[f(s, bigzones[i])])
    shs <- do_shifts(sub.rans[[1]], sub.rans[[2]], max.dist, min.dist)
    absolute.pairs <- lapply(shs, .absPos, sub.rans, sets)
    in.range <- lapply(absolute.pairs, .getInRange, sets, smallzones[i])

    return(in.range)
}

.joinZones <- function(xs, is, js)
{
    recSubset <- function(xs, i, j)
        xs[[i]][[j]]

    nestedLapply <- function(xs, ys, f)
        lapply(xs, function(x) lapply(ys, function(y) f(x, y)))

    joined <- nestedLapply(
        is, js,
        f <- function(i, j)
            unique(unlist(lapply(xs, recSubset, i, j)))
    )
    names(joined) <- is

    return(joined)
}

#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges start end
shifts <- function(setA, setB, win.size=10000, max.dist=74, min.dist=10)
{
    sets <- list(setA, setB)
    sets.range <- range(do.call(`c`, sets))

    sep <- floor(win.size) / 2

    zone.a <- .buildZoneRange(start(sets.range),
                              end(sets.range),
                              win.size)

    zone.b <- .buildZoneRange(start(sets.range) + sep,
                              end(sets.range),
                              win.size)

    if (any(sapply(list(zone.a, zone.b), is.null))) {
        return(list(left  = list(IRanges(), IRanges()),
                    right = list(IRanges(), IRanges()),
                    rest  = list(IRanges(), IRanges())))
    }

    a.shrinked <- .shrinker(zone.a, except.start=1)
    b.shrinked <- .shrinker(zone.b, except.end=length(zone.b))

    zone <- c(zone.a, zone.b)
    subzone <- c(a.shrinked, b.shrinked)

    by.zones <- lapply(seq_along(zone),
                       .doZone,
                       zone,
                       subzone,
                       sets,
                       max.dist,
                       min.dist)

    idxs <- .joinZones(by.zones, c("left", "right"), c(1, 2))

    left.shifts <- mapply(`[`, sets, idxs$left)
    right.shifts <- mapply(`[`, sets, idxs$right)

    f <- function (...)
        mapply(c, ..., SIMPLIFY=FALSE)
    rest <- mapply(`[`, sets, lapply(do.call(f, idxs), `-`))

    return(list(left  = left.shifts,
                right = right.shifts,
                rest  = rest))
}

.applyDistThresh <- function (rs, minDiff=10)
{
    ds <- lapply(rs, .dyadPos)
    diffs <- abs(do.call(`-`, ds))
    sel <- diffs >= minDiff

    pos.sel <- lapply(rs, `[`,  sel)
    neg.sel <- lapply(rs, `[`, !sel)

    list(properShifts=pos.sel, smallShifts=neg.sel)
}
