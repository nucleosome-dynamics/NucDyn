.doShifts <- function(sets, maxDiff)
{
    preciseDyadPos <- function(set)
        as.integer(((start(set) + end(set)) / 2) * 10)

    xDyad <- preciseDyadPos(sets[[1]])
    xSize <- as.integer(length(xDyad))
    yDyad <- preciseDyadPos(sets[[2]])
    ySize <- as.integer(length(yDyad))

    xLeft <- .initEmptyVect(xSize)
    yLeft <- .initEmptyVect(ySize)
    yRight <- yLeft
    xRight <- xLeft

    diff <- as.integer(maxDiff * 10)

    cOut <- .C("shifts",
               xDyad, xSize, yDyad, ySize,
               xLeft=xLeft, yLeft=yLeft, xRight=xRight, yRight=yRight,
               diff)

    return(list(left=list(cOut$xLeft, cOut$yLeft),
                right=list(cOut$xRight, cOut$yRight)))
}

.buildZoneRange <- function(start, end, width)
{
    npos <- ceiling((end - start) / width)
    zone.starts <- 0:(npos-1) * width + start
    zone.rans <- IRanges(start=zone.starts, width=width)
    return(zone.rans)
}

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

.getAbsolutePos <- function(idxs, xs, subran, wholeran)
{
    idxs <- as.integer(idxs)
    nidxs <- as.integer(length(idxs))

    xs <- as.integer(xs)
    xstart <- as.integer(start(subran))
    xend <- as.integer(end(subran))
    nxs <- as.integer(length(xs))

    whole.xstart <- as.integer(start(wholeran))
    whole.xend <- as.integer(end(wholeran))
    nwhole <- as.integer(length(wholeran))

    out <- as.integer(rep(0, nidxs))

    cOut <- .C("find_abs_pos",
               idxs, nidxs,
               xs, xstart, xend, nxs,
               whole.xstart, whole.xend, nwhole,
               out=out)

    return(cOut$out)
}

.absPos <- function(sh, sub.rans, sets)
{
    xs <- sh[[1]][as.logical(sh[[1]])]
    f <- function(a, b, c) .getAbsolutePos(xs, a, b, c)
    mapply(f, sh, sub.rans, sets, SIMPLIFY=FALSE)
}

.getInRange <- function(dyads, sets, ran)
{
    f <- function(a, xs)
        end(a) > start(xs) & start(a) < end(xs)

    pos <- do.call(`&`,
                   lapply(mapply(`[`,
                                 sets,
                                 dyads),
                          f,
                          ran))
    lapply(dyads, `[`, pos)
}

.doZone <- function(i, bigzones, smallzones, sets, dist)
{
    f <- function(a, xs)
        end(a) > start(xs) & start(a) < end(xs)

    sub.rans <- lapply(sets,
                       function(s) s[f(s, bigzones[i])])

    shs <- .doShifts(sub.rans, dist)

    absolute.pairs <- lapply(shs,
                             .absPos,
                             sub.rans,
                             sets)

    in.range <- lapply(absolute.pairs,
                       .getInRange,
                       sets,
                       smallzones[i])

    return(in.range)
}

.joinZones <- function(xs, is, js)
{
    recSubset <- function(xs, i, j)
        xs[[i]][[j]]

    nestedLapply <- function(xs, ys, f)
        lapply(xs,
               function(x) lapply(ys,
                                  function(y) f(x, y)))

    joined <- nestedLapply(
        is, js,
        f <- function(i, j)
            unique(unlist(lapply(xs, recSubset, i, j)))
    )
    names(joined) <- is

    return(joined)
}

.shifts <- function(sets, win.size=10000, dist=74)
{
    sets.range <- range(do.call(`c`, sets))

    sep <- floor(win.size) / 2

    zone.a <- .buildZoneRange(start(sets.range),
                              end(sets.range),
                              win.size)

    zone.b <- .buildZoneRange(start(sets.range) + sep,
                              end(sets.range),
                              win.size)

    a.shrinked <- .shrinker(zone.a, except.start=1)
    b.shrinked <- .shrinker(zone.b, except.end=length(zone.b))

    zone <- c(zone.a, zone.b)
    subzone <- c(a.shrinked, b.shrinked)

    by.zones <- lapply(seq_along(zone), .doZone, zone, subzone, sets, dist)

    idxs <- .joinZones(by.zones,
                       c("left", "right"),
                       c(1, 2))

    left.shifts <- mapply(`[`, sets, idxs$left)
    right.shifts <- mapply(`[`, sets, idxs$right)

    rest <- mapply(`[`,
                   sets,
                   lapply(do.call(function(...) mapply(c,
                                                       ...,
                                                       SIMPLIFY=FALSE),
                                  idxs),
                          `-`))

    return(list(left=left.shifts,
                right=right.shifts,
                rest=rest))
}
