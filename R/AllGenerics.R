#!/usr/bin/env Rscript

setGeneric(
    "nucleosomeDynamics",
    function(setA, setB, ...) standardGeneric("nucleosomeDynamics")
)

setGeneric(
    "findHotspots",
    function(dyn, ...) standardGeneric("findHotspots")
)

setGeneric(
    "plotDynamics",
    function(dyn, ...) standardGeneric("plotDynamics")
)

#setGeneric(
#    "applyThreshold",
#    function(hs, threshold, scale) standardGeneric("applyThreshold")
#)

setGeneric(
    "applyThreshold",
    function(hs, threshold, ...) standardGeneric("applyThreshold")
)
