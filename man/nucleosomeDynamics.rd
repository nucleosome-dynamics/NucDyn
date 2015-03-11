\name{nucleosomeDynamics}
\alias{nucleosomeDynamics}
\alias{nucleosomeDynamics,IRanges-method}
\alias{nucleosomeDynamics,GRanges-method}
\alias{nucleosomeDynamics,RangedData-method}
\title{
    Run nucleosomeDynamics to compare two sets.
}
\description{
    This is the main function in nucleosomeDynamics. It allows to compare the
    reads of two NGS experiments of nucleosome coverage.
}
\usage{
    \S4method{nucleosomeDynamics}{IRanges}(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140, maxDiff=74)
    \S4method{nucleosomeDynamics}{GRanges}(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140, maxDiff=74, mc.cores=1)
    \S4method{nucleosomeDynamics}{RangedData}(setA, setB, maxLen=170, equalSize=FALSE, roundPow=5, readSize=140, maxDiff=74, mc.cores=1)
}
\arguments{
    \item{setA}{
        Reads of the first experiment in \code{IRanges}, \code{RangedData},
        or \code{GRanges}. The format should the the same of the one in setB.
    }
    \item{setB}{
        Reads of the second experiment in \code{IRanges}, \code{RangedData},
        or \code{GRanges}. The format should the the same of the one in setA.
    }
    \item{maxLen}{
        Reads longer than this number will be filtered out.
    }
    \item{equalSize}{
        If set to \code{TRUE}, all sets will be set to the same length,
        conserving their original dyad position. Use it if the reads in your
        sets have differences in length (i.e. due to differences in the
        digestion) that you are not interested in.
    }
    \item{roundPow}{
        When equalSize is \code{FALSE}, the start and end of each read will be
        rounded to a power of this number to allow a more granular analysis.
        Set it to 0 if you don't want this granularization.
    }
    \item{readSize}{
        Length to which all reads will be set in case \code{equalSize} is
        \code{TRUE}.  It is ignored when \code{equalSize} is set to
        \code{FALSE}.
    }
    \item{maxDiff}{
        Maximum distance between the dyads of two reads that allows them to
        still be considered a "shift".
    }
    \item{mc.cores}{
        If \code{parallel} support, the number of cores available. This option
        is only used if the provided sets are from more than one chromosome.
    }
}
\value{
    An object of class \code{NucDyn} with two attributes: \code{set.a},
    \code{set.b}, each containing a \code{GRangesList} object with information
    about how reads where classified in each input dataset.
    The names of the \code{GRangesList} objects the following (with some
    exceptions*):
    \itemize{
        \item{originals}: All reads present for that set before the analysis.
        \item{coinciding}: Reads that have an exact match in the other set and
                           are considered to be the same read.
        \item{same.start}: Reads whose match in the other set starts at the
                           same position but ends in a different one. The are
                           considered to be a consequence of experimental
                           differences.
        \item{same.end}: Reads whose match in the other set ends at the same
                         position but starts in a different one. They are
                         considered to be a consequence of experimental
                         differences.
        \item{containedA}: Reads in the set B which are totally contained by
                           the length of their match in the set A. Or reads in
                           the set A that totally contain the length of their
                           match in the set B. They are considered to be a
                           consequence of experimental differences.
        \item{containedB}: Reads in the set A which are totally contained by
                           the length of their match in the set B. Or reads in
                           the set B that totally contain the length of their
                           match in the set A. They are considered to be a
                           consequence of experimental differences.
        \item{left.shifts}: Shifts that suffer an upstream shift when
                            considering transition from the set A to the set B.
        \item{right.shifts}: Shifts that suffer a downstream shift when
                             considering transition from the set A to the set
                             B.
        \item{indels}: Reads that are not present in the other set due to an
                       insertion or deletion event. If they are in set A, they
                       are considered a deletion and if they are in set B, they
                       are considered an insertion.
        \item{unpaired}: Reads that are uniformly and randomly removed from the
                         set with more reads to account for differences in the
                         number of reads between sets.
    }
    *: the types "same.start", "same.end", "containedA" and "containedB" will
       not be present if \code{equalSize} is set to \code{TRUE}.
}
\author{
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{manip}
