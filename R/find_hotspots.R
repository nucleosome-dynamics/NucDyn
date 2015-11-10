.subsetReads <- function (reads, range)
{  # Given a ranged data of reads and a range in the chormosome,
   # return the read subset in that range
    subset <- start(reads) >= range[[1]] & end(reads) <= range[[2]]
    subRange <- reads[subset]
    return(subRange)
}

.subsetCover <- function (reads, range)
{  # Given a ranged data of reads and a range in the chromosome,
   # return the coverage of the reads in that range
    width <- length(range[1]:range[2])
    subRange <- .subsetReads(reads, range)
    cov <- coverage(subRange, shift=(-range[1]), width=width)
    return(as.numeric(cov))
}

.splitBySign <- function (xs)
{
    # Convert a numeric vector into a list of two numeric vectors.
    # The first one will have the originaly positive numbers set to zero and
    # the second one, the originaly negative ones.
    a.prof <- xs
    a.prof[a.prof > 0] <- 0
    a.prof <- abs(a.prof)

    b.prof <- xs
    b.prof[b.prof < 0] <- 0

    list(a=a.prof, b=b.prof)
}

.myFilter <- function (xs)
    # Attempt to apply a FFT filter, return NULL if an error occurs
    tryCatch({
        filterFFT(xs, pcKeepComp=0.01, useOptim=TRUE)
    }, error=function(e) {
        NULL
    })

.setThreshold <- function (threshold="60%", xs)
{   # Handle a threshold possibly given as a string representing a precent
    if (is.null(xs)) {
        return(NULL)
    }
    if (is.character(threshold)) {
        t <- quantile(xs[xs != 0], as.numeric(sub("%", "", threshold))/100)
    } else {
        t <- threshold
    }
    if (t < 1) {
        return(1)
    } else {
        return(t)
    }
}

.sumWithNulls <- function (x, y)
{   # Sum two vectors and handle appropriately if one of them is NULL
    xcheck <- !is.null(x)
    ycheck <- !is.null(y)

    if (xcheck && ycheck) {
        x + y
    } else if (xcheck && !ycheck) {
        x
    } else if (!xcheck && ycheck) {
        y
    } else {
        NULL
    }
}

.handleNulls <- function(f)
    # Higer order function to properly handle NULL inputs
    function(x, ...)
        if (is.null(x) || (class(x) == "data.frame" && nrow(x) == 0)) {
            NULL
        } else {
            f(x, ...)
        }

.whichRan <- function(xs, r)
{   # Given a vector of positions and an IRange, return in which range each
    # position is found. Implemented in C.

    xs <- as.integer(xs)
    starts <- as.integer(start(r))
    ends <- as.integer(end(r))
    n <- as.integer(length(xs))
    m <- as.integer(length(r))
    out <- .initEmptyVect(n)

    cOut <- .C("which_ran", xs, n, starts, ends, m, out=out)$out

    #return(r[cOut])
    return(cOut)
}

.makeVectsEqual <- function(x, y)
{   # Make two numerical vectors the same size by appending zeros to the
    # shorter one
    vs <- list(x, y)
    lens <- sapply(vs, length)
    d <- max(lens) - min(lens)
    i <- which.min(lens)
    vs[[i]] <- c(vs[[i]], rep(0, d))
    return(vs)
}

.getHsRanges <- function(xs)
{   # Get all the hotspots ranges. Each range will be a non-zero range
    # separated by local minima
    nums <- c(xs, 0, 0)
    lims <- which(diff(sign(diff(nums))) > 0)
    all.ranges <- IRanges(start=c(1,
                                  lims[-length(lims)] + 2),
                          end=lims)
    return(all.ranges[nums[start(all.ranges)] != 0])
}

.buildDf <- .handleNulls(function(p, f)
    data.frame(coord=p, nreads=round(f[p])))

.applyThresh <- .handleNulls(function(df, t) {
    df[is.na(df$nreads), "nreads"] <- 0
    df[df$nreads > t, ]
})

.addRanges <- .handleNulls(function(df, cov) {
    set.range <- .getHsRanges(cov)
    idxs <- .whichRan(df$coord, set.range)
    df <- df[as.logical(idxs), ]
    rs <- set.range[idxs]
    df$start <- start(rs)
    df$end <- end(rs)
    df
})

.addTypes <- .handleNulls(function(df, n) {
    df$type <- n
    df
})

.peakCoupledParams <- function (ranA, ranB, range, threshold, names)
{   # Given two ranged data of reads that are coupled (ex: left and right
    # shifts, inserts and evictions), return the peaks and their scores

    rans <- lapply(list(ranA,
                        ranB),
                   .subsetReads,
                   range)

    covs <- do.call(.makeVectsEqual,
                    lapply(rans,
                           coverage))

    diff <- covs[[2]] - covs[[1]]
    by.sign <- .splitBySign(diff)
    filtered <- lapply(by.sign, .myFilter)

    peak.dfs <- mapply(
        function(f, n, t) {
            p <- .handleNulls(peakDetection)(f, threshold=0, score=FALSE)
            df <- .buildDf(p, f)
            df <- .applyThresh(df, t)
            df <- .addRanges(df, f)
            df <- .addTypes(df, n)
            df
        },
        filtered, as.list(names),
        SIMPLIFY=FALSE,
        MoreArgs=list(.setThreshold(threshold,
                        do.call(.sumWithNulls,
                                unname(filtered))))
    )

    df <- do.call(rbind, peak.dfs)

    if (is.null(df) || nrow(df) == 0) {
        return(data.frame(chrom  = numeric(0),
                          coord  = numeric(0),
                          type   = numeric(0),
                          nreads = numeric(0),
                          start  = numeric(0),
                          end    = numeric(0)))
    }

    ordered.df <- df[order(df$coord), ]
    return(ordered.df)
}

.quantileNormalize <- function(x, y)
{
    getReorderer <- function(x)
        order(c(1:length(x))[order(x)])
    ref <- sort((x + y) / 2)
    list(ref[getReorderer(x)],
         ref[getReorderer(y)])
}

.indelPeaks <- function(dyn, range, threshold="60%")
{
    reads <- lapply(dyn$originals,
                    function(x) .subsetReads(unique(.rmLongReads(x,
                                                                 174)),
                                             range=range))

    #covs <- lapply(dyn$indels, coverage)
    covs <- lapply(reads, coverage)

    norm.covs <- do.call(.quantileNormalize,
                         lapply(do.call(.makeVectsEqual,
                                        covs),
                                as.vector))

    filtered <- lapply(norm.covs,
                       filterFFT,
                       pcKeepComp=0.01,
                       useOptim=TRUE)

    diff <- do.call(`-`, filtered)

    insert.prof <- diff
    insert.prof[insert.prof > 0] <- 0
    insert.prof <- -insert.prof
    delete.prof <- diff
    delete.prof[delete.prof < 0] <- 0

    set.range <- .getHsRanges(abs(diff))

    thresh <- .setThreshold(threshold, abs(diff))

    ips <- peakDetection(insert.prof, threshold=thresh, score=FALSE)
    ips <- data.frame(peak=ips, score=insert.prof[ips])

    dps <- peakDetection(delete.prof, threshold=thresh, score=FALSE)
    dps <- data.frame(peak=dps, score=delete.prof[dps])

    idxs <- .whichRan(ips$peak, set.range)
    ips <- ips[as.logical(idxs), ]
    rs <- set.range[idxs]
    ips$start <- start(rs)
    ips$end <- end(rs)

    idxs <- .whichRan(dps$peak, set.range)
    dps <- dps[as.logical(idxs), ]
    rs <- set.range[idxs]
    dps$start <- start(rs)
    dps$end <- end(rs)

    ips$type <- "INCLUSION"
    dps$type <- "EVICTION"

    df <- rbind(ips, dps)
    sorted.df <- df[order(df$peak), ]

    return(with(sorted.df,
                data.frame(coord  = peak,
                           type   = type,
                           nreads = round(score),
                           start  = start,
                           end    = end)))
}

.typeNamer <- function(typeA, typeB, mergeKind)
{
    if (mergeKind == "shift-shift") {
        if (typeA == "SHIFT -" && typeB == "SHIFT +") {
            return("DISPERSION")
        } else if (typeA == "SHIFT +" && typeB == "SHIFT -") {
            return("CONCENTRATION")
        } else {
            return(NULL)
        }

    } else if (mergeKind == "shift-shift_again") {
        if (typeA == "OPENING -" && typeB == "OPENING +") {
            return("OPENING <>")
        } else if (typeA == "CLOSING +" && typeB == "CLOSING -") {
            return("CLOSING <>")
        } else {
            return(NULL)
        }

    } else if (mergeKind == "shift-indel") {
        forbidden <- ((typeA == "DISPERTION"    && typeB == "INCLUSION") ||
                      (typeA == "CONCENTRATION" && typeB == "EVICTION"))
        if (forbidden) {
            return(NULL)
        }

        if (typeA == "DISPERTION" || typeA == "CONCENTRATION") {
            direction <- "<>"
        } else if (typeA == "SHIFT +") {
            direction <- "+"
        } else if (typeB == "SHIFT -") {
            direction <- "-"
        } else {
            return(NULL)
        }

        if (typeB == "EVICTION") {
            kind <- "OPENING"
        } else if (typeB == "INCLUSION") {
            kind <- "CLOSING"
        } else {
            return(NULL)
        }

        return(paste(kind, direction, sep=" "))
    }
}

.merger <- function(overlap, peaksA, peaksB, mergeKind, same.magnitude,
                    nuc.peaks=NULL)
{   # Function that returns which new peak to create and which ones
    # to delete after a merge. Used in an lapply in mergeShifts and
    # in mergeShiftsIndels

    .isSameNuc <- function(coord1, coord2, peaks)
    {
        whichPeak <- function(pos, peaks)
            which(start(peaks) <= pos &
                  end(peaks) >= pos)

        nucs <- unlist(lapply(list(coord1,
                                   coord2),
                              whichPeak,
                              peaks))

        return(length(nucs) == 2 && length(unique(nucs)) == 1)
    }

    h1 <- queryHits(overlap)
    h2 <- subjectHits(overlap)
    e1 <- peaksA[h1, ]
    e2 <- peaksB[h2, ]

    # in shift mergings, only look at adjacent peaks belonging to the same
    # nucleosome
    if (mergeKind == "shift-shift") {
        not.adjacent <- h2 - h1 != 1  # TRUE if reads are not adjacent
        different.nuc <- !.isSameNuc(e1$coord, e2$coord, nuc.peaks)
        if (not.adjacent || different.nuc) {
            return(NULL)
        }
    }

    mergedPair <- rbind(e1, e2)

    # check if same magnitude
    compNReads <- mergedPair$nreads
    if (max(compNReads) / min(compNReads) <= same.magnitude) {
        newType <- .typeNamer(e1$type, e2$type, mergeKind)
        if (is.null(newType)) {
            return(NULL)
        }
        return(list(a2rm    = h1,
                    b2rm    = h2,
                    newPeak = with(mergedPair,
                                   data.frame(coord  = round(mean(coord)),
                                              type   = newType,
                                              nreads = sum(nreads),
                                              start  = min(start),
                                              end    = max(end)))))
    } else {
        return(NULL)
    }
}

.mergeShifts <- function (shift.peaks, nuc.width, mergeKind="shift-shift",
                          same.magnitude, nuc.peaks)
{   # Merge shifts
    if (nrow(shift.peaks) == 0) {
        return(shift.peaks)
    }

    shift.ran <- IRanges(start = shift.peaks$coord - round(nuc.width/2),
                         width = nuc.width)
    ovlp <- findOverlaps(shift.ran, ignoreSelf=TRUE, ignoreRedundant=TRUE)

    if (length(ovlp) > 0) {
        toRmAndAdd <- lapply(1:length(ovlp),
                             function(i) .merger(ovlp[i], shift.peaks,
                                                 shift.peaks, mergeKind,
                                                 same.magnitude,
                                                 nuc.peaks=nuc.peaks))
        shifts2rm <- unlist(lapply(toRmAndAdd,
                                   function(i) c(i$a2rm, i$b2rm)))
        newPeaks <- do.call("rbind",
                            lapply(toRmAndAdd,
                                   function(i) i$newPeak))

        if (length(shifts2rm) > 0) {
            shift.peaks <- shift.peaks[-shifts2rm, ]
            shift.peaks <- rbind(shift.peaks, newPeaks)
        }
    }
    return(shift.peaks)
}

.mergeShiftsIndels <- function (shift.peaks, indel.peaks, nuc.width,
                                mergeKind="shift-indel", same.magnitude)
{  # Merge shifts with indels
    if (nrow(shift.peaks) > 0 && nrow(indel.peaks) > 0) {
        shift.ran <- IRanges(start = shift.peaks$coord - round(nuc.width/2),
                             width = nuc.width)
        indel.ran <- IRanges(start = indel.peaks$coord - round(nuc.width/2),
                             width = nuc.width)
    } else {
        return(list(shifts=shift.peaks, indels=indel.peaks))
    }

    ovlp <- findOverlaps(shift.ran, indel.ran)

    if (length(ovlp) > 0) {
        toRmAndAdd <- lapply(1:length(ovlp),
                             function(i) .merger(ovlp[i],
                                                 shift.peaks,
                                                 indel.peaks,
                                                 mergeKind,
                                                 same.magnitude))
        shift2rm <- unlist(lapply(toRmAndAdd,
                                  function(i) i$a2rm))
        indel2rm <- unlist(lapply(toRmAndAdd,
                                  function(i) i$b2rm))
        newPeaks <- do.call("rbind",
                            lapply(toRmAndAdd,
                                   function(i) i$newPeak))
        if (length(shift2rm) > 0) {
            shift.peaks <- shift.peaks[-shift2rm, ]
        }
        shift.peaks <- rbind(shift.peaks, newPeaks)
        if (length(indel2rm) > 0) {
            indel.peaks <- indel.peaks[-indel2rm, ]
        }
    }
    return(list(shifts=shift.peaks, indels=indel.peaks))
}

.combinePeaks <- function(shift.peaks, indel.peaks, nuc.width, same.magnitude,
                          nuc.peaks)
{   # Merge peaks
    # find overlaps of shifts
    shift.peaks <- .mergeShifts(shift.peaks,
                                nuc.width=nuc.width,
                                same.magnitude=same.magnitude,
                                nuc.peaks=nuc.peaks)

    # find overlaps between indels and (new) shifts
    merged <- .mergeShiftsIndels(shift.peaks, indel.peaks, nuc.width=nuc.width,
                                 same.magnitude=same.magnitude)
    shift.peaks <- merged$shifts
    indel.peaks <- merged$indels

    # another round in the shifts for possible mid-distance combinations
    shift.peaks <- .mergeShifts(shift.peaks, mergeKind="shift-shift_again",
                                nuc.width=nuc.width,
                                same.magnitude=same.magnitude,
                                nuc.peaks=nuc.peaks)

    return(list(shifts=shift.peaks, indels=indel.peaks))
}

.getReadsInvolved <- function(coords, readsA, readsB)
{
    myFun <- "reads_involved"

    coords <- as.integer(coords)
    coordNum <- as.integer(length(coords))

    startsA <- as.integer(start(readsA))
    endsA <- as.integer(end(readsA))
    readNumA <- as.integer(length(startsA))

    startsB <- as.integer(start(readsB))
    endsB <- as.integer(end(readsB))
    readNumB <- as.integer(length(startsB))

    countOut <- as.double(replicate(coordNum, 0))

    .C(myFun, coords, coordNum,
       startsA, endsA, readNumA,
       startsB, endsB, readNumB,
       out=countOut)$out
}

.getPeaks <- function(reads, nuc.width)
{   # to do a nucleosome peak calling
    cov <- coverage.rpm(RangedData(reads))[[1]]
    filtered <- filterFFT(cov, pcKeepComp=0.01)
    peaks <- peakDetection(filtered, threshold="25%", score=FALSE,
                           width=nuc.width)
}

.findInRange <- function (dyn, range, nuc.width=120, combined=TRUE,
                          same.magnitude=2, threshold="60%")
{
    if (is.null(range)) {
        wholeRange <- range(do.call("c", dyn$originals))
        range <- sapply(c(start, end), function(f) f(wholeRange))
    }

    shift.peaks <- .peakCoupledParams(dyn$left.shifts[[1]],
                                      dyn$right.shifts[[1]],
                                      range,
                                      threshold,
                                      c("SHIFT -",
                                        "SHIFT +"))

    indel.peaks <- .indelPeaks(dyn, range, threshold=threshold)

    contained.peaks <- .peakCoupledParams(dyn$containedA[[1]],
                                          dyn$containedB[[2]],
                                          range,
                                          threshold,
                                          c("CONTAINED BinA",
                                            "CONTAINED AinB"))

    if (combined) {
        nuc.peaks <- .getPeaks(dyn$originals[[1]], nuc.width=nuc.width)
        merged <- .combinePeaks(shift.peaks,
                                indel.peaks,
                                nuc.width=nuc.width,
                                same.magnitude=same.magnitude,
                                nuc.peaks=nuc.peaks)
        shift.peaks <- merged$shifts
        indel.peaks <- merged$indels
    }

    all <- rbind(shift.peaks, indel.peaks, contained.peaks)
    all <- all[order(all$coord), ]

    # Relative number of reads involved in the hotspot
    numReads <- mean(unlist(lapply(c(1, 2),
                            function (n)
                                length(.subsetReads(dyn$originals[[n]],
                                                    range)))))
    #all$freads <- round((all$nreads/numReads) * 100) / 100
    all$freads <- all$nreads/numReads

    # Fraction of reads involved in position
    readsInvolved <- .getReadsInvolved(all$coord,
                                       dyn$originals[[1]],
                                       dyn$originals[[2]])
    all$hreads <- round((all$nreads/readsInvolved) * 100) / 100

    return(all)
}

setMethod(
    "findHotspots",
    signature(dyn="NucDyn"),
    function (dyn, range=NULL, chr=NULL, nuc.width=120, combined=TRUE,
              same.magnitude=2, threshold="60%", mc.cores=1) {
        setA <- set.a(dyn)
        setB <- set.b(dyn)

        chrs <- levels(seqnames(setA)$originals)

        chrIter <- function(chr, range=NULL) {
            message(paste("Starting", chr))
            f <- function(x) ranges(x[seqnames(x) == chr])
            chrDyn <- mapply(list, f(setA), f(setB), SIMPLIFY=FALSE)
            hs <- .findInRange(chrDyn,
                               range=range,
                               nuc.width=nuc.width,
                               combined=combined,
                               same.magnitude=same.magnitude,
                               threshold=threshold)
            if (nrow(hs)) {
                hs$chr <- chr
            }
            rownames(hs) <- NULL
            message(paste(chr, "done"))
            hs
        }

        if (is.null(chr)) {
            hsLs <- .xlapply(chrs, chrIter, mc.cores=mc.cores)
            return(do.call("rbind", hsLs))
        } else if (chr %in% chrs) {
            return(chrIter(chr, range=range))
        } else {
            stop("chromosome ", chr, " not found")
        }
    }
)
