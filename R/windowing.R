getSplitterIdxs <- function (n, wins, shift=0) {
    a <- seq(from=1+shift, to=n+1, by=wins)
    i <- a - floor(wins/2)
    j <- a + ceiling(wins/2) - 1
    i[i < 1] <- 1
    j[j > n] <- n
    list(i=i, j=j)
}

splitWithIdxs <- function (x, idxs)
    mapply(function(a, b) x[a:b],
           idxs[["i"]],
           idxs[["j"]],
           SIMPLIFY=FALSE)

windowSplitter <- function (x, wins=10000) {
    n <- length(x)
    unshifted.idxs <- getSplitterIdxs(n, wins, 0)
    shifted.idxs   <- getSplitterIdxs(n, wins, floor(wins/2))
    unshifted <- splitWithIdxs(x, unshifted.idxs)
    shifted   <- splitWithIdxs(x, shifted.idxs)
    list(unshifted=unshifted, shifted=shifted)
}

getShrinkerIdxs <- function (n, wins, shift=0) {
    init.idxs <- getSplitterIdxs(n, wins, shift)[["i"]]
    a <- seq(from=1+shift, to=n+1, by=wins)
    i <- a - ceiling(floor(wins/2)/2)
    j <- a + floor(ceiling(wins/2)/2) - 1
    i[i < 1] <- 1
    j[j > n] <- n
    i <- i - init.idxs + 1
    j <- j - init.idxs + 1
    list(i=i, j=j)
}

windowShrinker <- function (xs, n, wins, shift=0) {
    idxs <- getShrinkerIdxs(n, wins, shift)
    mapply(function (x, i, j) x[i:j],
           xs,
           idxs[["i"]],
           idxs[["j"]],
           SIMPLIFY=FALSE)
}

windowJoiner <- function (xs, ys, n, wins=10000) {
    xs <- windowShrinker(xs, n, wins)
    ys <- windowShrinker(ys, n, wins, floor(wins/2))

    if (length(xs) < length(ys)) {
        xs <- c(xs, list(c()))
    }
    if (length(xs) > length(ys)) {
        ys <- c(ys, list(c()))
    }

    unlist(mapply(c, xs, ys, SIMPLIFY=FALSE))
}

mapplyer <- function (f, x)
    do.call(function (...) mapply(f, ..., SIMPLIFY=FALSE), x)

doBySplitting <- function(f, wins, ...) {
    xs <- list(...)
    n <- min(sapply(xs, length))
    xs <- lapply(xs, `[`, 1:n)

    splitted <- lapply(xs, windowSplitter, wins=wins)

    aa <- lapply(splitted, `[[`, 1)
    bb <- lapply(splitted, `[[`, 2)
    raa <- mapplyer(f, aa)
    rbb <- mapplyer(f, bb)

    rejoined <- windowJoiner(raa, rbb, n, wins=wins)
    return(rejoined)
}
