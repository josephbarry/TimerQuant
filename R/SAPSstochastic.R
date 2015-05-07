simulatedRatio <- function(T1, T2, hLife, sigmaAdd, p=1, E=0) {
    FP1 <- x1fretFP1ss(p, log(2)/T1, log(2)/T2, log(2)/hLife, E=E)
    FP2 <- x1ss(p, log(2)/T2, log(2)/hLife)
    FP1 <- FP1+rnorm(1, 0, sigmaAdd)
    FP2 <- FP2+rnorm(1, 0, sigmaAdd)
    r <- FP2/FP1
    return(r)
}

simulatedSignal <- function(T1, T2, TA, TB, sigmaAdd, p=1, E=0)
    log2(simulatedRatio(T1, T2, TB, sigmaAdd, p, E)/
        simulatedRatio(T1, T2, TA, sigmaAdd, p, E))

simulatedSignalN <- function(T1, T2, TA, TB, sigmaAdd, N, p=1, E=0) {
    d <- sapply(seq_len(N), function(i) suppressWarnings(simulatedSignal(T1, T2,
        TA, TB, sigmaAdd, p, E)))
    return(d)
}

