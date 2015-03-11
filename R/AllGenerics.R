#!/usr/bin/env Rscript

setGeneric(
    "nucleosomeDynamics",
    function(setA, setB, equalSize=FALSE, mc.cores=1, ...)
        standardGeneric("nucleosomeDynamics")
)

setGeneric(
    "findHotspots",
    function(dyn, range=c(), chr=NULL, nuc.width=120, combined=TRUE,
             same.magnitude=2, threshold="60%", mc.cores=1, useOptim=FALSE)
        standardGeneric("findHotspots")
)
