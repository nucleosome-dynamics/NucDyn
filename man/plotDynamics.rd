\name{plotDynamics}
\alias{plotDynamics}
\alias{plotDynamics,NucDyn-method}
\title{
    Find hotspots in a NucDyn object.
}
\description{
    Find hotspots from a given nucleosome dynamics.
}
\usage{
    \S4method{plotDynamics}{NucDyn}(dyn, plot.range=NULL, chr=NULL, dyn.name="Dyn", expA.name="Ref 1", expB.name="Ref 2", norm.factor=1, ...)
}
\arguments{
    \item{dyn}{
        \code{NucDyn} object with the dynamic to plot.
    }
    \item{plot.range}{
        Range from the \code{NucDyn} object to plot. If not specified, the
        whole set will be plotted. If the dynamics contains more than one
        chromosome, they will appear concatenated in the plot.
    }
    \item{chr}{
        Chromosome from the \code{NucDyn} object to plot. If not specified, all
        chromosomes will appear plotted concatenated.
    }
    \item{dyn.name}{
        Name to be given to the dyanamics that will be displayed in the plot.
    }
    \item{expA.name}{
        Name to be given to the first data set of the dyanamics that will be
        displayed in the plot.
    }
    \item{expB.name}{
        Name to be given to the second data set of the dyanamics that will be
        displayed in the plot.
    }
    \item{norm.factor}{
        Normalization factor between ref1 and ref2. Use it to visualize both
        coverages profiles on a similar scale if one of them has a
        significantly higher coverage.
    }
    \item{\dots}{
        Other parameters passed to \code{\link{plot}} function.
    }
}
\details{
    Plotting offers a visual representation of the dynamics between two
    reference states.

    Coverage profile for ref1 will be shown as a solid grey background and
    coverage profile for ref2 will be shown as a dotted profile.
    Superimposed, arrows showing read shifts (upstream in blue and downstream
    in red).
    Indels will appear as a small coverage profile at the bottom part of the
    plot (insertions in green and deletions in red).
}
\value{
    Void
}
\author{
    Oscar Flores \email{oflores@mmb.pcb.ub.es}\cr
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{hplot}
