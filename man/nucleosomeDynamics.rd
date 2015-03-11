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
        or \code{GRanges}. The format should the same of the one in setB.
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
        sets have differences in length (ie, due to differences in the
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
\details{
    The aim of NucleosomeDynamics is to infer "movement" (with direction and
    magnitude) of the reads between two reference nucleosome maps. In contrast
    with a simple coverage difference, NucleosomeDynamics can tell how the
    reads change between two different experiments. This is useful to analyze
    regions where fine regulatory role of the nucleosomes is suspected to
    happen.

    This method is based on the idea that reads in a reference state (ref1)
    should match those in another reference state (ref2) after applying a few
    shifts and/or indels. Both ref1 and ref2 need to be experimental nucleosome
    maps, either from the same sample with different conditions or from
    different samples.

    Before the analysis read positions are granularized by setting their starts
    and ends to numbers that are a power of \code{roundPow} (in the case where
    \code{equalSize=FALSE}) or set to the length of \code{readSize} (if
    \code{equalSize=TRUE}).

    Then, we look for a match to each read in ref1 in another read of ref2
    using a specific deffinition for what a "match" is. After all possible
    matches have been found, we set those reads apart and we look for matches
    in the remaining reads using a different definition for what a "match" is.
    Then, to account for the possibility that one sample has a higher coverage
    than the other, randomly picked reads are removed from the dataset with
    more reads to that they are the same size.
    After this, all the remaining reads are considered indels.

    Since they are tried sequentially, the definition of each type of match
    implies that the previous definitions tried are do not hold.
    The different types of matches, in the order in which they are tried are:
    \itemize{
        \item{Coinciding}{
            Reads that start and end in the exact same position in both sets.
        }
        \item{Same start}{
            Reads that start in the same position in both sets but end at a
            different one.
        }
        \item{Same end}{
            Reads that end in the same position in both sets but start at a
            different one.
        }
        \item{Contained}{
            Reads from one set that are contained or contain reads from the
            other set. For a read to be contained by another, it has to start
            at a more upstream position but end in a more downstream position
            than the second read.
        }
        \item{Shifts}{
            Reads whose dyads are at a maximum distance of \code{maxDist}
            (default is 74 bp.).
        }
    }
    If \code{equalSize=TRUE}, all reads are forced to the same size, to the
    match types "Same start", "Same end" and "Contained" do not apply.

    In an attempt to find the optimum pairing for the shifts, we use dynamic
    programming approach. It is done in such a way that the score is inversely
    proportional to the dyad distance (to favor shortest possible distance),
    but with a very high gap penalty (to favor more pairs at longer distance
    rather than less closer pairs) and a penalty of \code{-Inf} for distances
    higher than \code{maxDist} so that those don't happen at all.
}
\value{
    An object of class \code{NucDyn} with two attributes: \code{set.a},
    \code{set.b}, each containing a \code{GRangesList} object with information
    about how reads where classified in each input dataset.
    The names of the \code{GRangesList} objects the following (with some
    exceptions*):
    \itemize{
        \item{originals}{
            All reads present for that set before the analysis.
        }
        \item{coinciding}{
            Reads that have an exact match in the other set and are considered
            to be the same read.
        }
        \item{same.start}{
            Reads whose match in the other set starts at the same position but
            ends in a different one. The are considered to be a consequence of
            experimental differences.
        }
        \item{same.end}{
            Reads whose match in the other set ends at the same position but
            starts in a different one. They are considered to be a consequence
            of experimental differences.
        }
        \item{containedA}{
            Reads in the set B which are totally contained by the length of
            their match in the set A. Or reads in the set A that totally
            contain the length of their match in the set B. They are considered
            to be a consequence of experimental differences.
        }
        \item{containedB}{
            Reads in the set A which are totally contained by the length of
            their match in the set B. Or reads in the set B that totally
            contain the length of their match in the set A. They are considered
            to be a consequence of experimental differences.
        }
        \item{left.shifts}{
            Shifts that suffer an upstream shift when considering transition
            from the set A to the set B.
        }
        \item{right.shifts}{
            Shifts that suffer a downstream shift when considering transition
            from the set A to the set B.
        }
        \item{indels}{
            Reads that are not present in the other set due to an insertion or
            deletion event. If they are in set A, they are considered a
            deletion and if they are in set B, they are considered an
            insertion.
        }
        \item{unpaired}{
            Reads that are uniformly and randomly removed from the set with
            more reads to account for differences in the number of reads
            between sets.
        }
    }
    *: the types "same.start", "same.end", "containedA" and "containedB" will
       not be present if \code{equalSize} is set to \code{TRUE}.
}
\author{
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{manip}
