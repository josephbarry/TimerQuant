\name{plotPrimordiumProfile}
\alias{plotPrimordiumProfile}
\title{
Plot Primordium Profile
}
\description{
Visualizes primordium signal as a function of position with median and median absolute deviation across samples.
}
\usage{
plotPrimordiumProfile(x, add, ylab, lwd, cex.lab, cex.axis, xlim, ylim, main, col, lty, alpha)
}
\arguments{
  \item{x}{
A matrix of data where rows are samples and columns are sequential positions.
}
  \item{add}{
A logical indicating whether or not to add to the existing plot.
}
  \item{ylab}{
The y-axis label.
}
  \item{lwd}{
Integer specifying width of lines.
}
  \item{cex.lab}{
Integer specifying size of labels.
}
  \item{cex.axis}{
Integer specifying size of axis labels.
}
  \item{xlim}{
An optional vector of length 2 specifying the limits for the x-axis.
}
  \item{ylim}{
An optional vector of length 2 specifying the limits for the y-axis.
}
  \item{main}{
Plot title.
}
  \item{col}{
Line colour.
}
  \item{lty}{
Style of line
}
  \item{alpha}{
A numeric between zero and one specifying the level of transparency for the shaded region.
}
}
\value{
Produces a plot of signal vs position summarizing across multiple primordium samples.
}
\author{
Joseph D. Barry
}
\keyword{plotPrimordiumProfile}
