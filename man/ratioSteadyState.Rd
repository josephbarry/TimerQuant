\name{ratioSteadyState}
\alias{ratioSteadyState}
\title{
analytic function ratioSteadyState
}
\description{
Steady-state analytic solution to one-step model for the ratio of mature to 
non-mature fluorophores.
}
\usage{
ratioSteadyState(T1, T2, halfLife, E=0, f=1)
}
\arguments{
\item{T1}{
Maturation time of fluorescent protein 1 (FP1, fast maturing).
}
\item{T2}{
Maturation time of fluorescent protein 2 (FP2, slow maturing).
}
\item{halfLife}{
Protein half-life.
}
\item{E}{
FRET value representing transfer from FP1 to FP2.
}
\item{f}{
f=f2/f1, the ratio of prefactors relating the number of molecules to 
fluorescence intensity for each fluorescence channel.
}
}
\value{
A numeric specifying the model steady-state solution for the given parameters.
}
\author{
Joseph D. Barry
}
\examples{
halfLifeSeq <- seq(1, 2000, by=0.1)
plot(halfLifeSeq, ratioSteadyState(T1=5, T2=100, halfLife=halfLifeSeq), 
    type="l", lwd=2, ylim=c(0, 1), xlab="tFT half-life (min)", 
    ylab="Steady-state ratio", cex.lab=1.4, log="x", col="red")
}
\keyword{ratioSteadyState}
