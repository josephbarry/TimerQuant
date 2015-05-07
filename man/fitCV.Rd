\name{fitCV}
\alias{fitCV}
\title{
fitCV
}
\description{
Fits a smoothing line to coefficient of variation profiles.
}
\usage{
fitCV(x, scaleLog10)
}
\arguments{
\item{x}{
A dataframe with columns \code{T} (FP maturation time) and \code{CV}
(coefficient of variation of timer signal).
}
\item{scaleLog10}{
A logical indicating whether the points are spaced on the log10 scale or not.
}
}
\value{
A dataframe containing fitted values and the minimum CV of the profile.
}
\examples{
vignette(topic="genPaperFigures", package="TimerQuant")
}
\author{
Joseph D. Barry
}
\keyword{fitCV}
