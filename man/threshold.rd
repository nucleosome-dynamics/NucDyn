\name{threshold}
\alias{applyThreshold}
\alias{applyThreshold,data.frame,ThresholdByType-method}
\alias{applyThreshold,data.frame,VariableThreshold-method}
\alias{applyThreshold,data.frame,character-method}
\alias{applyThreshold,data.frame,integer-method}
\alias{applyThreshold,data.frame,numeric-method}
\alias{getVariableThreshold}
\alias{variableThreshFromHss}
\title{
    Deal with threshold objects for hotspots.
}
\description{
    Deal with threshold objects for hotspots.
}
\usage{
    \S4method{applyThreshold}{data.frame,numeric}(hs, threshold)
    \S4method{applyThreshold}{data.frame,integer}(hs, threshold)
    \S4method{applyThreshold}{data.frame,character}(hs, threshold)
    \S4method{applyThreshold}{data.frame,VariableThreshold}(hs, threshold, scale=NULL)
    \S4method{applyThreshold}{data.frame,ThresholdByType}(hs, threshold)
}
\arguments{
    \item{hs}{
        Hotspots returned by findHotspots.
    }
    \item{threshold}{
        Threshold object.
    }
    \item{scale}{
        Scaling on a variable threshold.
    }
}
\author{
    Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{manip}
