#!/usr/bin/env Rscript

setGeneric(
    "nucleosomeDynamics",
    function(setA, setB, ...) standardGeneric("nucleosomeDynamics")
)

setGeneric(
    "findHotspots",
    function(dyn, ...) standardGeneric("findHotspots")
)
