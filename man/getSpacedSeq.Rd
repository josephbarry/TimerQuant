\name{getSpacedSeq}
\alias{getSpacedSeq}
\title{
Get Spaced Sequence
}
\description{
Return points nicely spaced for on the log10 scale.
}
\usage{
getSpacedSeq(x, n)
}
\arguments{
\item{x}{
A vector of two numbers containing the minimum and maximum of the desired 
sequence.
}
\item{n}{
The desired length of the sequence to be returned.
}
}
\value{
A sequence of numbers with appropriate spacing for the log10 scale.
}
\examples{
getSpacedSeq(c(1, 1000), n=10)
}
\author{
Joseph D. Barry
}
\keyword{getSpacedSeq}
