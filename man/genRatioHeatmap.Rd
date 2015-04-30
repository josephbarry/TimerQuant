\name{genRatioHeatmap}
\alias{genRatioHeatmap}
\title{
generate ratio heatmap
}
\description{
Visualize timer ratios as a function of either FP1 or FP2 maturation time, and protein half-life.
}
\usage{
genRatioHeatmap(tRangeFP, Tfixed, TA, TB, channel, E, f=1, n, ramp)
}
\arguments{
  \item{tRangeFP}{
Vector containing two numerics specifying the range of maturation times to display for the chosen fluorescence channel.
}
  \item{Tfixed}{
The maturation time of the FP that will remain fixed.
}
  \item{TA}{
tFT half-life in location A.
}
  \item{TB}{
tFT half-life in location B.
}
  \item{channel}{
Integer specifying fluorescence channel to be varied (1 or 2).
}
 \item{E}{
FRET value representing transfer from FP1 to FP2.
}
  \item{f}{
f=f2/f1, the ratio of prefactors relating the number of molecules to the fluorescence intensity.
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
\keyword{genRatioHeatmap}
