#!/usr/bin/env Rscript

VariableThreshold <- setClass("VariableThreshold",
                              representation(x="list",
                                             scale="numeric"))

setMethod(
    "show",
    signature = "VariableThreshold",
    definition = function (object) {
        cat("variable threshold\n")
        invisible(NULL)
    }
)

ThesholdByType <- setClass("ThresholdByType",
                           representation(shifts="numeric",
                                          indels="numeric",
                                          contained="numeric"))

setMethod(
    "show",
    signature = "ThresholdByType",
    definition = function (object) {
        cat("threshold by type\n")
        cat("shifts: ", object@shifts, "\n")
        cat("indels: ", object@indels, "\n")
        cat("contained: ", object@contained, "\n")
    }
)
