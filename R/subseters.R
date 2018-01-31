.separateGroups <- function(sets, subsetList)
{
    subA <- find_subs(subsetList[[1]])
    subB <- find_subs(subsetList[[2]])

    seleA <- sets[[1]][subA]
    seleB <- sets[[2]][subB]

    nonSeleA <- sets[[1]][!as.logical(subsetList[[1]])]
    nonSeleB <- sets[[2]][!as.logical(subsetList[[2]])]

    return(list(matches = list(seleA,    seleB),
                rest    = list(nonSeleA, nonSeleB)))
}

.separateShifts <- function(sets, subsetList)
{
    subLeftA <- find_subs(subsetList$left[[1]])
    subLeftB <- find_subs(subsetList$left[[2]])
    subRightA <- find_subs(subsetList$right[[1]])
    subRightB <- find_subs(subsetList$right[[2]])

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
