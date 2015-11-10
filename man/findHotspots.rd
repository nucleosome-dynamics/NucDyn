\name{findHotspots}
\alias{findHotspots}
\alias{findHotspots,NucDyn-method}
\title{
    Find hotspots in a NucDyn object.
}
\description{
    Find hotspots from a given nucleosome dynamics.
}
\usage{
    \S4method{findHotspots}{NucDyn}(dyn, range=NULL, chr=NULL, nuc.width=120, combined=TRUE, same.magnitude=2, threshold="60\%", mc.cores=1)
}
\arguments{
    \item{dyn}{
        \code{NucDyn} object with the dynamic to analyze.
    }
    \item{range}{
        Range from the \code{NucDyn} object to analyze. If not specified, the
        whole set will be analyzed.
    }
    \item{chr}{
        Chromosome from the \code{NucDyn} object to analyze. If not specified,
        all chromosomes will be analyzed.
    }
    \item{nuc.width}{
        Nucleosome width considered.
    }
    \item{combined}{
        If \code{TRUE}, nearby hotspots will be combined (see details for
        possible combinations). Otherwise, all the information about detected
        shifts and changes in coverage will be returnes without further
        processing.
    }
    \item{same.magnitude}{
        Only used if \code{combined=TRUE}. When combining two hotspots, this
        value is the maximum ratio between two spots to be considered with the
        same magnitude. Two hot spots with the same magnitude can be combined,
        but a large hotspot will not be merged with a much smaller one. See
        details.
    }
    \item{threshold}{
        \code{character} value giving a relative coverage percentile or
        \code{numeric} value giving an absolute threshold. Hotspots with less
        than this number of reads will be discarded.
    }
    \item{mc.cores}{
        If \code{parallel} support, the number of cores available. This option
        is only used if the provided sets are from more than one chromosome.
    }
}
\details{
    This function is aimed to help in the analysis of the nucleosome dynamics
    by pointing out these regions with relevant changes in the individual
    position of the nucleosomes.

    There are 4 types of basic hotspots, plus their combinations. Basic ones
    are:

    \itemize{
        \item{"SHIFT +"}: Translational movement of nucleosomes, downstream
                          (+).
        \item{"SHIFT -"}: Translational movement of nucleosomes, upstream (-).
        \item{"EVICTION"}: Nucleosome reads removed from that locus.
        \item{"INCLUSION"}: Nucleosome reads added to that locus.
    }

    If \code{combined=TRUE}, adjacent and overlapped hotspots will be combined
    as follows:

    \itemize{
        \item{"CONCENTRATION"}: "SHIFT +" followed by "SHIFT -".
        \item{"DISPERSION"}: "SHIFT -" followed by "SHIFT +".
        \item{"OPENING <>"}: "EVICTION" overlapped with "DISPERSION".
        \item{"OPENING +"}: "EVICTION" overlapped with "SHIFT +".
        \item{"OPENING -"}: "EVICTION" overlapped with "SHIFT -".
        \item{"CLOSING <>"}: "INCLUSION" overlapped with "CONCENTRATION".
        \item{"CLOSING +"}: "INCLUSION" overlapped with "SHIFT +".
        \item{"CLOSING -"}: "INCLUSION" overlapped with "SHIFT -".
    }

    As Translational and coverage changes can happen anywhere, only those
    involving a certain number of reads are reported. This number can by
    adjusted by the \code{threshold} parameter. If \code{threshold} is a
    \code{character} vector representing a percentage value (ie,
    \code{"60\%"}), this will be automatically converted to the absolute value
    given by the corresponding percentile of the coverage in the window. If,
    instead, \code{threshold} is a \code{numeric} value, this value will be
    used as absolute threhold.

    It two adjacent hotspots with shifts in opposite directions are detected
    but one of them is relatively small in comparison with the other, but will
    be reported as shifts, disregarding the value of \code{combined}. We
    consider two hotspots of the same magnitude if the ratio between the number
    of reads in one and the other is smaller than \code{same.magnitude}.
    This ratio is always performed by using the larger number as numerator and
    the smaller as denominator; therefore, \code{same.magnitude} must always be
    greater of equal than 1.

    For example, with \code{same.magnitude=2}, we consider that 25 reads
    shifting downstream followed bby 17 reads shifting upstream will be of the
    same magnitude (25/17 == 1.47 < 2) and we will annotate it as a
    "DISPERSION". In another example, if we have 25 shifts downstream followed
    by only 5 shifts upstream (25/5 == 5 > 2), both hotspot will be annotated
    as "SHIFT".
}
\value{
    A \code{data.frame} with the following columns:

    \itemize{
        \item{chrom:}{
            Chromosome name.
        }
        \item{coord:}{
            Genomic coordinates (average dyad position of affected
            nucleosomes).
        }
        \item{type:}{
            The type of the hotspot (as listed above).
        }
        \item{nreads:}{
            Number of reads involved in the hotspot.
        }
    }
}
\author{
    Oscar Flores \email{oflores@mmb.pcb.ub.es}\cr
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{manip}
