#!/usr/bin/env Rscript

setGeneric(
    "nucleosomeDynamics",
    function(setA, setB, equalSize=FALSE, mc.cores=1, ...) standardGeneric("nucleosomeDynamics")
)
