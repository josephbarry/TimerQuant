getSpacedSeq <- function(x, n) 10^seq(log10(x[1]), log10(x[2]), length=n)

getBreaks10 <- function(x) {
    lo <- floor(log10(min(x, na.rm=TRUE)))
    hi <- ceiling(log10(max(x, na.rm=TRUE)))
    breaks <- as.vector(10^(lo:hi) %o% c(1, 5))
    breaks <- breaks[order(breaks, decreasing=FALSE)]
    return(breaks)
}

solveForHeatPlot <- function(mTimeMin, mTimeMax, pAgeMin, pAgeMax,
    seqLength=110, channel=2, Tfixed=5, E=0, f=1) {
    proteinRange <- getSpacedSeq(c(log(2)/pAgeMin, log(2)/pAgeMax), seqLength)
    mRange <- getSpacedSeq(c(log(2)/mTimeMin, log(2)/mTimeMax), seqLength)
    gridSolve <- expand.grid(mRange, proteinRange)
    colnames(gridSolve) <- c("m", "protein")
    gridSolveList <- as.list(data.frame(t(gridSolve)))
    if (channel == 2) ratioSolve=sapply(gridSolveList,
        function(x) ratioSteadyState(T1=Tfixed, T2=log(2)/x[1], 
        halfLife=log(2)/x[2], E=E, f=f))
    else if (channel == 1) ratioSolve <- sapply(gridSolveList,
        function(x) ratioSteadyState(T2=Tfixed, T1=log(2)/x[1], 
        halfLife=log(2)/x[2], E=E, f=f))
    else stop("channel expected to be 1 or 2")
    gridSolve$ratio <- ratioSolve
    gridSolve$mTime <- log(2)/gridSolve$m
    gridSolve$proteinHlife <- log(2)/gridSolve$protein
    return(gridSolve)
}

genRatioHeatmap <- function(tRangeFP, Tfixed, TA, TB, channel=2, E=0, f=1, 
    n=110, ramp=rainbow(100)) {
    x <- solveForHeatPlot(mTimeMin=tRangeFP[1], mTimeMax=tRangeFP[2], 
        pAgeMin=TA,pAgeMax=TB, channel=channel, Tfixed=Tfixed, E=E, f=f, 
        seqLength=n)
    hmap <- qplot(x=x$mTime, y=x$proteinHlife, colour=x$ratio, shape=I(15),
        xlab=paste0("FP", channel, " maturation time (minutes)"),
        ylab=expression("protein half-life"~T[1/2]~"(minutes)"))+
        coord_cartesian(xlim=tRangeFP, ylim=c(TA, TB))+
        scale_x_log10(breaks=getBreaks10(tRangeFP))+
        scale_y_log10(breaks=getBreaks10(c(TA, TB)))+
        scale_linetype_discrete(name="a")+
        scale_colour_gradientn(colours=ramp)+
        labs(colour=expression("R"))+
        theme_bw()+
        theme(axis.title.x=element_text(size=32), 
            axis.title.y=element_text(size=32, vjust=1.5), 
            axis.text.x=element_text(size=24), 
            axis.text.y=element_text(size=24), 
            legend.key.size=unit(64, "points"), 
            legend.title=element_text(size=24),
            legend.text=element_text(size=24), panel.border=element_blank())
    return(hmap)
}

genTimeSteadyStateHeatmap <- function(tRangeFP2, tRangeHlife, n=110, 
    ramp=rainbow(100)) {
    FP2seq <- getSpacedSeq(tRangeFP2, n)
    hlifeseq <- getSpacedSeq(tRangeHlife, n)
    grid <- expand.grid(FP2seq, hlifeseq)
    colnames(grid) <- c("mTime", "hLife")
    gridList <- as.list(data.frame(t(grid)))
    grid$tss <- sapply(gridList, function(x) tss(log(2)/x[1], log(2)/x[2]))
    hmap <- qplot(x=grid$mTime, y=grid$hLife, colour=grid$tss, shape=I(15),
        xlab="FP2 maturation time (minutes)", 
        ylab="protein half-life (minutes)")+
        coord_cartesian(xlim=tRangeFP2, ylim=tRangeHlife)+
        scale_x_log10(breaks=getBreaks10(tRangeFP2))+
        scale_y_log10(breaks=getBreaks10(tRangeHlife))+
        scale_colour_gradientn(colours=ramp)+
        labs(colour=expression(T[ss]))+
        theme_bw()+
        theme(axis.title.x=element_text(size=32), 
            axis.title.y=element_text(size=32, vjust=1.5), 
            axis.text.x=element_text(size=24),
            axis.text.y=element_text(size=24),
            legend.key.size=unit(64, "points"),
            legend.title=element_text(size=24),
            legend.text=element_text(size=24), panel.border=element_blank())
    return(hmap)
}

ratioSteadyState <- function (T1, T2, halfLife, E=0, f=1) {
    k <- log(2)/halfLife
    m1 <- log(2)/T1
    m2 <- log(2)/T2
    return(f*(m1+k)/(m2*(1-E)+k)*m2/m1)
}

runChooseFP2App <- function() 
    runApp(system.file("extdata/chooseFP2App", package="TimerQuant"))

runTimerModellingApp <- function()
    runApp(system.file("extdata/timerModellingApp", package="TimerQuant"))

ratioTimeDependent <- function(T1, T2, halfLife, t, E=0, f=1) {
    k <- log(2)/halfLife
    m1 <- log(2)/T1
    m2 <- log(2)/T2
    ratio <- (exp(t*m1)*(-1+exp(k*t))*f*(k+m1)*(k+exp(t*m2)*
        (-k+(-1+exp(k*t))*m2)))/((-k+exp(t*m1)*(k-(-1+exp(k*t))*m1))*
        (E*k+exp(t*m2)*(-(-1+E+exp(k*t))*k+(-1+E)*(-1+exp(k*t))*m2)))
    return(ratio)
}

x0ss <- function(p, m, k) p/(k+m)

x1 <- function (p, m, k, t, f=1) 
    f*p/(k*(k+m))*exp(-(k+m)*t)*(k-exp(m*t)*(k-m*exp(k*t)+m))

x1ss <- function(p, m, k, f=1) f*p*m/(k*(k+m))

x1fretFP1 <- function(p, m1, m2, k, t, E=0, f=1)
    (exp(-t*(k+m1+m2))*f*p*(k+exp(t*m1)*(-k+(-1+exp(k*t))*m1))*((exp(t*m2)*
    (-1+exp(k*t))+E*(-1+exp(t*m2)))*k-(-1+E)*exp(t*m2)*(-1+exp(k*t))*m2))/
    ((-1+exp(k*t))*k*(k+m1)*(k+m2))

x1fretFP1ss <- function(p, m1, m2, k, E=0, f=1)
    f*p*m1*(k-(-1+E)*m2)/(k*(k+m1)*(k+m2))

tss <- function(m, k) 1/k+1/(k+m)+1/m*log(1+m/k)

inflectionSlope <- function(p, m, k, a) -(((-1+a)*m*((k+m)/k)^(-(k/m))*p)/
    (k+m))

inflectionOffset <- function(p, m, k, a)
    ((k+m)/k)^(1-k/m)*(m*(k+m)*((k+m)/k)^(k/m)*p+(-1+a)*p*(m*(2*k+m)+k*(k+m)*
    log((k+m)/k)))/(k+m)^3

signal <- function(T1, T2, TA, TB, E=0) {
    # note that production rate and brighness prefactors cancel out by taking 
    # the ratio so set to one for simplicity
    FP1A <- x1fretFP1ss(p=1, m1=log(2)/T1, m2=log(2)/T2, k=log(2)/TA, 
        E=E, f=1)
    FP2A <- x1ss(p=1, m=log(2)/T2, k=log(2)/TA, f=1)
    FP1B <- x1fretFP1ss(p=1, m1=log(2)/T1, m2=log(2)/T2, k=log(2)/TB, 
        E=E, f=1)
    FP2B <- x1ss(p=1, m=log(2)/T2, k=log(2)/TB, f=1)
    ra <- FP2A/FP1A
    rb <- FP2B/FP1B
    d <- log2(rb/ra)
    return(d)
}

fitCV <- function(x, scaleLog10=TRUE) {
    expectedNames <- c('Time', 'CV')
    if (! all(expectedNames %in% colnames(x))) 
        stop("expecting columns 'Time' and 'CV'")
    tSeq <- getSpacedSeq(range(x$Time, na.rm=TRUE), n=100)
    if (scaleLog10) x$Time <- log10(x$Time)
    poly <- lp(x$Time, deg=2)
    fit <- locfit(CV ~ Time, data=x)
    f <- function(x, scaleLog10) {
        if (scaleLog10) return(log10(x)) 
        else return(x)
    }
    pred <- data.frame(Time=tSeq, CV=predict(fit, 
        newdata=data.frame(Time=f(tSeq, scaleLog10))))
    fn <- approxfun(pred$Time, pred$CV)
    o <- optimize(fn, interval=range(pred$Time), maximum=FALSE)
    pred$FP2optimumTime <- o$minimum
    pred$FP2optimum <- o$objective
    return(pred)
}
