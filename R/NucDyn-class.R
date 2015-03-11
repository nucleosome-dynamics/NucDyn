#!/usr/bin/env Rscript

NucDyn <- setClass("NucDyn",
                   representation(set.a = "GRangesList",
                                  set.b = "GRangesList"))

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

setGeneric("set.a", function(object, ...) standardGeneric("set.a"))
setMethod("set.a", "NucDyn", function(object) object@set.a)
setGeneric("set.b", function(object, ...) standardGeneric("set.b"))
setMethod("set.b", "NucDyn", function(object) object@set.b)

