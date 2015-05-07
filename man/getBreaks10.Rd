\name{getBreaks10}
\alias{getBreaks10}
\title{
Get log10 breaks
}
\description{
Return breaks for each half-decade on the log10 scale, e.g. 1, 5, 10, 50, ...
}
\usage{
getBreaks10(x)
}
\arguments{
\item{x}{
A vector of numbers. Breaks will be calculated across the range of \code{x}.
}
}
\value{
A sequence of breaks useful for ticks or labels on the log10 scale.
}
\examples{
getBreaks10(c(1, 100))
}
\author{
Joseph D. Barry
}
\keyword{getBreaks10}
