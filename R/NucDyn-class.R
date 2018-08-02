#' Nucleosome Dynamics calculation.
#'
#' This class represents the result of running `nucleosomeDynamics` on two
#' reference sets.
#'
#' Intances of this class contain two `GRangesList` objects, each with
#' information about how reads where classified in each input dataset.  The
#' names of the `GRangesList` objects the following (with some exceptions*):
#'
#' * originals: All reads present for that set before the analysis.
#' * coinciding: Reads that have an exact match in the other set and are
#'   considered to be the same read.
#' * same.start: Reads whose match in the other set starts at the same position
#'   but ends in a different one. The are considered to be a consequence of
#'   experimental differences.
#' * same.end: Reads whose match in the other set ends at the same position but
#'   starts in a different one. They are considered to be a consequence of
#'   experimental differences.
#' * containedA: Reads in the set B which are totally contained by the length
#'   of their match in the set A. Or reads in the set A that totally contain
#'   the length of their match in the set B. They are considered to be a
#'   consequence of experimental differences.
#' * containedB: Reads in the set A which are totally contained by the length
#'   of their match in the set B. Or reads in the set B that totally contain
#'   the length of their match in the set A. They are considered to be a
#'   consequence of experimental differences.
#' * left.shifts: Shifts that suffer an upstream shift when considering
#'   transition from the set A to the set B.
#' * right.shifts: Shifts that suffer a downstream shift when considering
#'   transition from the set A to the set B.
#' * indels: Reads that are not present in the other set due to an insertion or
#'   deletion event. If they are in set A, they are considered a deletion and
#'   if they are in set B, they are considered an insertion.
#' * unpaired: Reads that are uniformly and randomly removed from the set with
#'   more reads to account for differences in the number of reads between sets.
#'
#' *: the types "same.start", "same.end", "containedA" and "containedB" will
#' not be present if in the call of `nucleosomeDynamics`, `equalSize=TRUE`.
#'
#' @return An instance of the class `NucDyn` representing the output of
#'   [nucleosomeDynamics()].
#'
#' @slot set.a `GRangesList` containing the `nucleosomeDynamics` of ref1
#' @slot set.b `GRangesList` containing the `nucleosomeDynamics` of ref2
#'
#' @author Ricard Illa \email{ricard.illa@@irbbarcelona.org}
#' @seealso [nucleosomeDynamics()]
#' @keywords manip
#' @importClassesFrom GenomicRanges GRangesList
#' @import doParallel Rcpp BiocGenerics
#'
NucDyn <- setClass(
    "NucDyn",
    representation(
        set.a = "GRangesList",
        set.b = "GRangesList"
    )
)

setMethod(
    "show",
    signature = "NucDyn",
    definition = function (object) {

        getLsName <- function (xs)
            paste(lapply(
                names(object@set.a),
                function(i)
                    paste0(i, ": ", length(xs[[i]]), " reads")
            ), collapse=" | ")

        cat("Set a:", "\n", sep="")
        cat(getLsName(object@set.a), "\n", sep="")
        cat("Set b:", "\n", sep="")
        cat(getLsName(object@set.b), "\n", sep="")

        invisible(NULL)
    }
)

#' Accessors to a NucDyn object
#'
#' @rdname nucdyn-accessors
#' @param x a NucDyn object
#' @return
#'   * {set.a} `GRangesList` with the reads paired from ref1
#'   * {set.b} `GRangesList` with the reads paired from ref2
#'
setGeneric("set.a", function (x) standardGeneric("set.a"))

#' @rdname nucdyn-accessors
setGeneric("set.b", function (x) standardGeneric("set.b"))

setMethod("set.a", "NucDyn", function (x) x@set.a)
setMethod("set.b", "NucDyn", function (x) x@set.b)
