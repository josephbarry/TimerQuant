\name{genTimeSteadyStateHeatmap}
\alias{genTimeSteadyStateHeatmap}
\title{
Time to reach steady-state
}
\description{
Visualize the time to reach steady-state as a function of FP2 maturation time and protein half-life.
}
\usage{
genTimeSteadyStateHeatmap(tRangeFP2, tRangeHlife, n, ramp)
}
\arguments{
  \item{tRangeFP2}{
Vector containing two numerics specifying the range of FP2 maturation times.
}
  \item{tRangeHlife}{
Vector containing two numerics specifying the range of protein half-lives.
}
  \item{n}{
Integer specifying the number of data points. Choose a higher \code{n} for a higher pixel density.
}
  \item{ramp}{
    Colour ramp, see \code{colorRampPalette} for more details.
}
}
\value{
Returns a ggplot2 heatmap.
}
\author{
Joseph D. Barry
}
\keyword{genTimeSteadyStateHeatmap}
