\name{signal}
\alias{signal}
\title{
Signal
}
\description{
Computes timer signal (without additive noise) for a set of model parameters.
}
\usage{
signal(T1, T2, TA, TB, E=0)
}
\arguments{
\item{T1}{
Maturation time of fluorescent protein 1 (fast maturing).
}
\item{T2}{
Maturation time of fluorescent protein 2 (slow maturing).
}
\item{TA}{
Minimum protein half-life.
}
\item{TB}{
Maximum protein half-life.
}
\item{E}{
FRET value representing transfer from FP1 to FP2.
}
}
\value{
A numeric specifying the timer signal.
}
\examples{
signal(T1=5, T2=60, TA=30, TB=180, E=0)
signal(T1=5, T2=60, TA=30, TB=180, E=0.5)
}
\author{
Joseph D. Barry
}
\keyword{signal}
