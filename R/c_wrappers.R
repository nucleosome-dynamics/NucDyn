.equals <- function(sets)
{
    myFun <- "equals"
    startA <- as.integer(start(sets[[1]]))
    endA <- as.integer(end(sets[[1]]))
    lenA <- as.integer(length(sets[[1]]))

    startB <- as.integer(start(sets[[2]]))
    endB <- as.integer(end(sets[[2]]))
    lenB <- as.integer(length(sets[[2]]))

    outA <- as.integer(rep(0, lenA))
    outB <- as.integer(rep(0, lenB))

    cOut <- .C(myFun,
               startA, endA, lenA,
               startB, endB, lenB,
               outA=outA, outB=outB)

    return(list(cOut$outA, cOut$outB))
}

.equalsAtDist <- function(sets, maxDist=5)
{
    myFun <- "same_at_dist"

    startsA <- as.integer(start(sets[[1]]))
    lenA <- as.integer(length(startsA))
    outA <- as.integer(rep(0, lenA))

    startsB <- as.integer(start(sets[[2]]))
    lenB <- as.integer(length(startsB))
    outB <- as.integer(rep(0, lenB))

    cOut <- .C(myFun, as.integer(maxDist),
               startsA, lenA,
               startsB, lenB,
               outA=outA, outB=outB)

    return(list(cOut$outA, cOut$outB))
}

.sameStart <- function(sets)
{
    myFun <- "same_start"
    startA <- as.integer(start(sets[[1]]))
    lenA <- as.integer(length(sets[[1]]))

    startB <- as.integer(start(sets[[2]]))
    lenB <- as.integer(length(sets[[2]]))

    outA <- as.integer(rep(0, lenA))
    outB <- as.integer(rep(0, lenB))

    cOut <- .C(myFun,
               startA, lenA,
               startB, lenB,
               outA=outA, outB=outB)

    return(list(cOut$outA, cOut$outB))
}

.sameEnd <- function(sets)
{
    myFun <- "same_end"
    endA <- as.integer(end(sets[[1]]))
    lenA <- as.integer(length(sets[[1]]))

    endB <- as.integer(end(sets[[2]]))
    lenB <- as.integer(length(sets[[2]]))

    outA <- as.integer(rep(0, lenA))
    outB <- as.integer(rep(0, lenB))

    cOut <- .C(myFun,
               endA, lenA,
               endB, lenB,
               outA=outA, outB=outB)

    return(list(cOut$outA, cOut$outB))
}

.contained <- function(sets)
{
    myFun <- "contained"
    startA <- as.integer(start(sets[[1]]))
    endA <- as.integer(end(sets[[1]]))
    lenA <- as.integer(length(sets[[1]]))

    startB <- as.integer(start(sets[[2]]))
    endB <- as.integer(end(sets[[2]]))
    lenB <- as.integer(length(sets[[2]]))

    outA <- as.integer(rep(0, lenA))
    outB <- as.integer(rep(0, lenB))

    cOut <- .C(myFun,
               startA, endA, lenA,
               startB, endB, lenB,
               outA=outA, outB=outB)

    return(list(cOut$outA, cOut$outB))
}

.preShifts <- function(sets, maxDiff)
{
    myFun <- "pre_shifts"
    startA <- as.integer(start(sets[[1]]))
    endA <- as.integer(end(sets[[1]]))
    lenA <- as.integer(length(sets[[1]]))

    startB <- as.integer(start(sets[[2]]))
    endB <- as.integer(end(sets[[2]]))
    lenB <- as.integer(length(sets[[2]]))

    leftA <- as.integer(rep(0, lenA))
    leftB <- as.integer(rep(0, lenB))

    rightA <- as.integer(leftA)
    rightB <- as.integer(leftB)

    cOut <- .C(myFun,
               startA, endA, lenA,
               startB, endB, lenB,
               leftA=leftA, leftB=leftB,
               rightA=rightA, rightB=rightB,
               as.integer(maxDiff))

    return(list(left=list(cOut$leftA, cOut$leftB),
                right=list(cOut$rightA, cOut$rightB)))
}

.shifts <- function (sets, shiftRan, maxDiff)
{
    myFun <- "shifts"
    setA <- sets[[1]]
    setB <- sets[[2]]

    xStart <- as.integer(start(setA))
    xEnd <- as.integer(end(setA))
    xSize <- as.integer(length(setA))

    yStart <- as.integer(start(setB))
    yEnd <- as.integer(end(setB))
    ySize <- as.integer(length(setB))

    ranStart <- as.integer(start(shiftRan))
    ranEnd <- as.integer(end(shiftRan))
    ranSize <- as.integer(length(shiftRan))

    xLeft <- as.integer(rep(0, xSize))
    yLeft <- as.integer(rep(0, ySize))
    xRight <- as.integer(rep(0, xSize))
    yRight <- as.integer(rep(0, ySize))

    cOut <- .C(myFun,
               xStart, xEnd, xSize,
               yStart, yEnd, ySize,
               ranStart, ranEnd, ranSize,
               xLeft=xLeft, yLeft=yLeft,
               xRight=xRight, yRight=yRight,
               as.double(maxDiff))

    return (list(left=list(cOut$xLeft, cOut$yLeft),
                 right=list(cOut$xRight, cOut$yRight)))
}

