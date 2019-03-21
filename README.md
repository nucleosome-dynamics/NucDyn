# NucleosomeDynamics

R package aimed at comparing the reads of two MNase-seq experiments for nucleosome positioning and detecting significant inclusions, evictions and shifts. 


## Installation
---------------

1- NucDyn depends on several R packages. To install them run:

&nbsp;

    source("https://bioconductor.org/biocLite.R")
    biocLite('nucleR')
    biocLite('dplyr')
    biocLite('IRanges')
    biocLite('GenomicRanges')
    biocLite('ShortRead')
    biocLite('doParallel')
    biocLite('ggplot2')
    biocLite('magrittr')


2- Download NucDyn from GitLab. 

`NucDyn_0.99.0.tar.gz`

The package has been developed and tested in R-3.5 (recommended version).

3- Install the package in R:
&nbsp;

   install.packages("NucDyn_0.99.0.tar.gz", repos = NULL)



## Usage
---------------

The main functionalities of NucDyn are performed with the functions  

* nucleosomeDynamics: finds differences between the two experiments at the fragment level
* findHotspots: combines fragment level information into hotspots of significant nucleosome changes 

An example of usage with the data provided with the package is explained below. For detailed explanation of the data used, how to interpret the results or how to run the functions with other data see the package vignette. 

1- Load the package in R 
&nbsp;

    library(NucDyn)

2- Load the example data for the two experimental conditions (see help for details)
&nbsp;
   
    data(readsG2_chrII)
    data(readsM_chrII)

3- Find differenecs at the fragment level
&nbsp;

    dyn <- nucleosomeDynamics(setA=readsG2_chrII, setB=readsM_chrII)
   
4- Combine fragment level results to detect hotspots of nucleosome movements
&nbsp;

    findHotspots(dyn, nuc_chrII)





## Developers
-------------

Diana Buitrago

Ricard Illa 

Diego Gallego
&nbsp;

&nbsp;

*Molecular Modeling and Bioinformatics Group.*


