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

.peakCoupledParams <- function (ranA, ranB, range, threshold, trim.ends,
                                useOptim)
{  # Given two ranged data of reads that are coupled (ex: left and right
   # shifts, inserts and evictions), return the peaks
    width <- length(range[1]:range[2])

    a.cov <- .subsetCover(ranA, range)
    b.cov <- .subsetCover(ranB, range)

    filtered <- tryCatch({
        filterFFT(b.cov - a.cov, useOptim=useOptim, pcKeepComp=0.01)
    }, error=function(e) {
        NULL
    })

    if (is.null(filtered)) {
        return(NULL)
    }

    a.filtered <- filtered
    a.filtered[a.filtered > 0] <- 0
    a.filtered <- abs(a.filtered)

    b.filtered <- filtered
    b.filtered[b.filtered < 0] <- 0

    if (is.character(threshold)) {
        percent <- as.numeric(sub("%", "", threshold))
        threshold <- quantile(abs(filtered), percent/100)
    }

    a.peaks <- peakDetection(a.filtered, threshold=threshold, score=FALSE)
    b.peaks <- peakDetection(b.filtered, threshold=threshold, score=FALSE)

    myPeaks <- c(a.peaks, b.peaks)

    if (!is.null(myPeaks) & length(myPeaks) > 0) {
        # join, order and add scores to peaks
        ordered <- order(myPeaks)
        myPeaks <- data.frame(peaks=myPeaks[ordered],
                              score=filtered[myPeaks][ordered])

        # only peaks with at least 1 whole read
        myPeaks[is.na(myPeaks$score), "score"] <- 0
        myPeaks <- myPeaks[abs(myPeaks$score) > 1, ]

        # remove peaks at the ends and keep only peks with at least
        #1 whole read
        subset <- !(myPeaks$peaks < trim.ends |
                    myPeaks$peaks > (width)-trim.ends) &
                  abs(myPeaks$score) > 1
        myPeaks <- myPeaks[subset, ]

        return(myPeaks)
    } else {
        return(NULL)
    }
}

.getSense <- function(scoreVect, nameVect)
{  # Use a vector of scores to return the appropriate names according to
   # whether the score is positive or negative
    return(unlist(lapply(
        scoreVect,
        function(score) {
            if (sign(score) == +1) {
                return(nameVect[1])
            } else {
                return(nameVect[2])
            }
        }
    )))
}

.shiftPeaks <- function(dyn, range, threshold="60%", trim.ends=25, useOptim)
{  # Return a data frame for the peaks of the shifts with its
   # corresponding classification name
    myPeaks <- .peakCoupledParams(dyn$left.shifts[[1]], dyn$right.shifts[[1]],
                                  range, threshold, trim.ends,
                                  useOptim=useOptim)
    if (!is.null(myPeaks) && nrow(myPeaks)) {
        sense <- .getSense(myPeaks$score, c("+", "-"))

        return(data.frame(coord  = myPeaks$peak+range[1],
                          type   = paste("SHIFT", sense),
                          nreads = round(abs(myPeaks$score))))
    } else {
        return(data.frame(chrom=numeric(0), coord=numeric(0),
                          type=numeric(0), nreads=numeric(0)))
    }
}

.containedPeaks <- function(dyn, range, threshold="60%", trim.ends=25,
                            useOptim)
{
    myPeaks <- .peakCoupledParams(dyn$containedA[[1]], dyn$containedB[[2]],
                                  range, threshold, trim.ends,
                                  useOptim=useOptim)
    if (!is.null(myPeaks) && nrow(myPeaks)) {
        sense <- .getSense(myPeaks$score, c("AinB", "BinA"))

        return(data.frame(coord  = myPeaks$peak+range[1],
                          type   = paste("CONTAINED", sense),
                          nreads = round(abs(myPeaks$score))))
    } else {
        return(data.frame(chrom=numeric(0), coord=numeric(0),
                          type=numeric(0), nreads=numeric(0)))
    }
}

.indelPeaks <- function(dyn, range, threshold="60%", trim.ends=25, useOptim)
{  # Return a data frame for the peaks of the indels with its
   # corresponding classification name
    myPeaks <- .peakCoupledParams(dyn$indels[[1]], dyn$indels[[2]],
                                  range, threshold, trim.ends,
                                  useOptim=useOptim)
    if (!is.null(myPeaks) && nrow(myPeaks)) {
        sense <- .getSense(myPeaks$score, c("INCLUSION", "EVICTION"))

        return(data.frame(coord  = myPeaks$peak+range[1],
                          type   = sense,
                          nreads = round(abs(myPeaks$score))))
    } else {
        return(data.frame(chrom=numeric(0), coord=numeric(0),
                          type=numeric(0), nreads=numeric(0)))
    }
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
        whichPeak <- function(peaks, pos)
            which(start(peaks) <= pos &
                  end(peaks) >= pos)

        nuc1 <- whichPeak(peaks, coord1)
        nuc2 <- whichPeak(peaks, coord2)

        return(length(nuc1) == 1 &&
               length(nuc2) == 1 &&
               nuc1 == nuc2)
    }

    h1 <- queryHits(overlap)
    h2 <- subjectHits(overlap)
    e1 <- peaksA[h1, ]
    e2 <- peaksB[h2, ]

    # in shift mergings, only look at adjacent peaks belonging to the same
    # nucleosome
    if (mergeKind == "shift-shift" &&
          (h2 - h1 != 1 ||  # TRUE if reads are not adjacent
           !.isSameNuc(e1$coord, e2$coord, nuc.peaks))) {
        return(NULL)
    }

    mergedPair <- rbind(e1, e2)

    # check if same magnitude
    compNReads <- mergedPair$nreads
    if (max(compNReads) / min(compNReads) <= same.magnitude) {
        newType <- .typeNamer(e1$type, e2$type, mergeKind)
        if (is.null(newType)) {
            return(NULL)
        }
        newPeak <- data.frame(coord  = round(mean(mergedPair$coord)),
                              type   = newType,
                              nreads = sum(mergedPair$nreads))
        return(list(a2rm=h1, b2rm=h2, newPeak=newPeak))
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
                             function(i) .merger(ovlp[i], shift.peaks,
                                                 indel.peaks, mergeKind,
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
    shift.peaks <- .mergeShifts(shift.peaks, nuc.width=nuc.width,
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
                          same.magnitude=2, threshold="60%", useOptim=FALSE)
{
    if (is.null(range)) {
        wholeRange <- range(do.call("c", dyn$originals))
        range <- sapply(c(start, end), function(f) f(wholeRange))
    }

    shift.peaks <- .shiftPeaks(dyn, range, threshold=threshold,
                               useOptim=useOptim)
    indel.peaks <- .indelPeaks(dyn, range, threshold=threshold,
                               useOptim=useOptim)
    contained.peaks <- .containedPeaks(dyn, range, threshold=threshold,
                                       useOptim=useOptim)

    if (combined) {
        nuc.peaks <- .getPeaks(dyn$originals[[1]], nuc.width=nuc.width)
        merged <- .combinePeaks(shift.peaks, indel.peaks, nuc.width=nuc.width,
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
    function (dyn, range=c(), chr=NULL, nuc.width=120, combined=TRUE,
              same.magnitude=2, threshold="60%", useOptim=FALSE, mc.cores=1) {
        setA <- set.a(dyn)
        setB <- set.b(dyn)

        chrs <- levels(seqnames(setA)$originals)

        chrIter <- function(chr) {
            message(paste("Starting", chr))

            f <- function(x) ranges(x[seqnames(x) == chr])
            chrDyn <- mapply(list, f(setA), f(setB), SIMPLIFY=FALSE)

            hs <- .findInRange(chrDyn, range=c(), nuc.width=nuc.width,
                               combined=combined,
                               same.magnitude=same.magnitude,
                               threshold=threshold, useOptim=useOptim)
            if (nrow(hs)) {
                hs$chr <- chr
            }
            message(paste(chr, "done"))
            hs
        }

        hsLs <- .xlapply(chrs, chrIter, mc.cores=mc.cores)

        return(do.call("rbind", hsLs))
    }
)
