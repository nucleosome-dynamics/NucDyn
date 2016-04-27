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
