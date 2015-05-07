\name{ratioTimeDependent}
\alias{ratioTimeDependent}
\title{
analytic function ratioTimeDependent
}
\description{
Steady-state analytic solution to one-step model for the ratio of mature to non-mature fluorophores.
}
\usage{
ratioTimeDependent(T1, T2, halfLife, t, E=0, f=1)
}
\arguments{
\item{T1}{
Maturation time of fluorescent protein 1 (FP1, fast maturing).
}
\item{T2}{
Maturation time of fluorescent protein 2 (FP2, slow maturing)
}
\item{halfLife}{
Protein half-life.
}
\item{t}{
Time, which must be non-negative.
}
\item{E}{
FRET value representing energy transfer from FP1 to FP2.
}
\item{f}{
f=f2/f1, the ratio of prefactors relating the number of molecules to fluorescence intensity.
}
}
\value{
A numeric specifying the model time-dependent solution for the given parameters.
}
\author{
Joseph D. Barry
}
\examples{
tSeq <- seq(0.1, 300, by=0.1)
plot(tSeq, ratioTimeDependent(T1=5, T2=100, halfLife=30, t=tSeq, E=0, f=1), type="l", lwd=2, 
    xlab="time (min)", ylab="ratio", cex.lab=1.4, col="black", ylim=c(0, 0.3))
points(tSeq, ratioTimeDependent(T1=5, T2=100, halfLife=30, t=tSeq, E=0.4, f=1), type="l", lwd=2, 
    col="red")
abline(h=ratioSteadyState(T1=5, T2=100, halfLife=30, E=0, f=1), lty=2, col="black")
abline(h=ratioSteadyState(T1=5, T2=100, halfLife=30, E=0.4, f=1), lty=2, col="red")
}
\keyword{ratioTimeDependent}
