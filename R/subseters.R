.findSubs <- function(subseter)
{
    myFun <- "find_subs"
    if (!is.loaded(myFun)) {
        dyn.load("find_subs.so")
    }

    subset <- as.integer(subseter)
    idxArray <- as.integer(sort(subseter[as.logical(subseter)]))

    inLen <- as.integer(length(subset))
    outLen <- as.integer(length(idxArray))

    out <- as.integer(rep(0, outLen))

    cOut <- .C(myFun, subset, idxArray, inLen, outLen, out=out)

    return(cOut$out)
}

.separateGroups <- function(sets, subsetList)
{
    subA <- .findSubs(subsetList[[1]])
    subB <- .findSubs(subsetList[[2]])

    seleA <- sets[[1]][subA]
    seleB <- sets[[2]][subB]

    nonSeleA <- sets[[1]][!as.logical(subsetList[[1]])]
    nonSeleB <- sets[[2]][!as.logical(subsetList[[2]])]

    return(list(matches = list(seleA,    seleB),
                rest    = list(nonSeleA, nonSeleB)))
}

.separateShifts <- function(sets, subsetList)
{
    subLeftA <- .findSubs(subsetList$left[[1]])
    subLeftB <- .findSubs(subsetList$left[[2]])
    subRightA <- .findSubs(subsetList$right[[1]])
    subRightB <- .findSubs(subsetList$right[[2]])

    seleLeftA <- sets[[1]][subLeftA, ]
    seleLeftB <- sets[[2]][subLeftB, ]
    seleRightA <- sets[[1]][subRightA, ]
    seleRightB <- sets[[2]][subRightB, ]

    nonSubA <- !(as.logical(subsetList$left[[1]]) |
                 as.logical(subsetList$right[[1]]))
    nonSubB <- !(as.logical(subsetList$left[[2]]) |
                 as.logical(subsetList$right[[2]]))

    nonSeleA <- sets[[1]][nonSubA, ]
    nonSeleB <- sets[[2]][nonSubB, ]

    return(list(left  = list(seleLeftA,  seleLeftB),
                right = list(seleRightA, seleRightB),
                rest  = list(nonSeleA,   nonSeleB)))
}
