\name{NucDyn-class}
\docType{class}

\alias{class:NucDyn}
\alias{NucDyn-class}

% accessors
\alias{set.a}
\alias{set.b}
\alias{set.a,NucDyn-method}
\alias{set.b,NucDyn-method}

\title{
    Nucleosome Dynamics calculation.
}
\description{
    This class represents the result of running \code{nucleosomeDynamics} on
    two reference sets.
}
\details{
    Intances of this class contain two \code{GRangesList} objects, each  with
    information about how reads where classified in each input dataset.
    The names of the \code{GRangesList} objects the following (with some
    exceptions*):
    \itemize{
        \item{originals:}{
            All reads present for that set before the analysis.
        }
        \item{coinciding:}{
            Reads that have an exact match in the other set and are considered
            to be the same read.
        }
        \item{same.start:}{
            Reads whose match in the other set starts at the same position but
            ends in a different one. The are considered to be a consequence of
            experimental differences.
        }
        \item{same.end:}{
            Reads whose match in the other set ends at the same position but
            starts in a different one. They are considered to be a consequence
            of experimental differences.
        }
        \item{containedA:}{
            Reads in the set B which are totally contained by the length of
            their match in the set A. Or reads in the set A that totally
            contain the length of their match in the set B. They are considered
            to be a consequence of experimental differences.
        }
        \item{containedB:}{
            Reads in the set A which are totally contained by the length of
            their match in the set B. Or reads in the set B that totally
            contain the length of their match in the set A. They are considered
            to be a consequence of experimental differences.
        }
        \item{left.shifts:}{
            Shifts that suffer an upstream shift when considering transition
            from the set A to the set B.
        }
        \item{right.shifts:}{
            Shifts that suffer a downstream shift when considering transition
            from the set A to the set B.
        }
        \item{indels:}{
            Reads that are not present in the other set due to an insertion or
            deletion event. If they are in set A, they are considered a
            deletion and if they are in set B, they are considered an
            insertion.
        }
        \item{unpaired:}{
            Reads that are uniformly and randomly removed from the set with
            more reads to account for differences in the number of reads
            between sets.
        }
    }
    *: the types "same.start", "same.end", "containedA" and "containedB" will
       not be present if in the call of \code{nucleosomeDynamics},
       \code{equalSize=TRUE}.
}
\value{
    An instance of the class \code{NucDyn} with the following slots:

    \item{set.a}{
        \code{GRangesList} containing the \code{nucleosomeDynamics} of ref1.
    }
    \item{set.b}{
        \code{GRangesList} containing the \code{nucleosomeDynamics} of ref2.
    }
}
\section{Accessor methods}{
    \describe{
        \item{}{
            \code{set.a(x)}: The dataset 1 from \code{x}.
        }
        \item{}{
            \code{set.b(x)}: The dataset 2 from \code{x}.
        }
    }
}
\author{
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\seealso{
    \code{\link{nucleosomeDynamics}}
}
\keyword{manip}
