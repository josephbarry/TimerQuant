\name{analyticSolutions}
\alias{tss}
\alias{x0ss}
\alias{x1}
\alias{x1ss}
\alias{x1fretFP1}
\alias{x1fretFP1ss}
\title{
Analytic model solutions for fluorescence intensity
}
\description{
Time-dependent and steady-state analytic solution to one-step model for number 
of mature fluorophores. Where \code{f} is given as a parameter the returned 
value is transformed from a molecular population into a fluorescence intensity.
For the function names, 0 refers to the dark population of non-mature 
fluorophores, and 1 to the mature, fluorescent population. 'ss' indicates 
steady-state solutions. \code{tss} is the time required to reach steady-state.
}
\usage{
tss(m, k)
x0ss(p, m, k)
x1(p, m, k, t, f=1)
x1ss(p, m, k, f=1)
x1fretFP1(p, m1, m2, k, t, E=0, f=1)
x1fretFP1ss(p, m1, m2, k, E=0, f=1)
}
\arguments{
\item{p}{
Protein production rate (molecules produced per unit time).
}
\item{m}{
Maturation rate of fluorophore, which can be for either FP1 or FP2 (convert to 
maturation time with log(2)/m).
}
\item{m1}{
Maturation rate of FP1.
}
\item{m2}{
Maturation rate of FP2.
}
\item{k}{
Protein degradation rate (convert to half-life with log(2)/k).
}
\item{t}{
Time (must be non-negative).
}
\item{E}{
FRET coefficient representing energy transfer from FP1 to FP2.
}
\item{f}{
Proportionality factor relating intensity to the number of molecules. When equal 
to one then the readout is number of molecules directly.
}
}
\value{
A numeric specifying the model solution for the given parameters.
}
\author{
Joseph D. Barry
}
\examples{
t0 <- seq(0.001, 1000, by=0.1)
plot(t0, x1(p=10, m=log(2)/5, k=log(2)/100, t=t0), type="l", col="darkgreen", 
    lwd=2, xlab="Time (min)", ylab="Number of mature fluorophores", cex.lab=1.4)
}
\keyword{x1}
