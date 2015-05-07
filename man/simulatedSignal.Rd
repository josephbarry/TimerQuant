\name{simulatedSignal}
\alias{simulatedSignal}
\alias{simulatedSignalN}
\alias{simulatedRatio}
\title{
Simulated Timer Signal
}
\description{
Additive error model for timer signal.
}
\usage{
simulatedRatio(T1, T2, hLife, sigmaAdd, p, E)
simulatedSignal(T1, T2, TA, TB, sigmaAdd, p, E)
simulatedSignalN(T1, T2, TA, TB, sigmaAdd, N, p, E)
}
\arguments{
\item{T1}{
Maturation time of fluorescent protein 1 (FP1, fast maturing).
}
\item{T2}{
Maturation time of fluorescent protein 2 (FP2, slow maturing).
}
\item{hLife}{
Protein half-life.
}
\item{TA}{
Minimum protein half-life.
}
\item{TB}{
Maximum protein half-life.
}
\item{sigmaAdd}{
Standard deviation of normal distribition from which noise terms are drawn.
}
\item{p}{
Protein production rate (molecules produced per unit time).
}
\item{E}{
FRET value representing energy transfer from FP1 to FP2.
}
\item{N}{
Number of simulation realizations.
}
}
\value{
Returns simulated values for ratios or timer signal.
}
\examples{
vignette(topic="genPaperFigures", package="TimerQuant")
}
\author{
Joseph D. Barry
}
\keyword{simulatedSignal}
