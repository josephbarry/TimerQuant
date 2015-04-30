library(TimerQuant)
library(grid)
library(ggplot2)
library(reshape2)
library(shiny)
library(shinyBS)
library(dplyr)
library(deSolve)
n <- 1000

myTheme <- theme_bw()+theme(title=element_text(size=18), text=element_text(size=12), 
    axis.title.y=element_text(vjust=1.5), plot.title=element_text(size=12))

modelFun <- function(model, t0, t1, tmax, r0, dr) {
    fn <- switch(model,
        constant=function(t) rep(r0, length(t)),
        linear=function(t) ifelse(t < t0, r0, ifelse(t > t1, r0+dr, r0+(t-t0)*dr/(t1-t0))),
        step=function(t) stepfun(c(t0, t1), c(r0, r0+dr, r0))(t))
    return(fn)
}

solveModel <- function(x01, x02, tseq, m1, m2, p0, k0, pfn, kfn) {
    eqs <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), list(c(dX1=p(t)-m*X1-k(t)*X1, 
            dX2=m*X1-k(t)*X2)))
    }
    fp1 <- ode(y=x01, times=tseq, func=eqs, parms=c(p=pfn, m=m1, k=kfn))
    fp2 <- ode(y=x02, times=tseq, func=eqs, parms=c(p=pfn, m=m2, k=kfn))
    df <- data.frame(time=fp1[, 1], FP1=fp1[, 3], FP2=fp2[, 3])
    df$Ratio <- df$FP2/df$FP1
    return(df)
}

shinyServer(function(input, output, session) {
    
    output$modelPlot <- renderPlot({
        tseq <- seq(0, input$Tf, length=n)
        pfn <- modelFun(model=input$prodModel, t0=input$startProd, t1=input$endProd, 
            tmax=input$Tf, r0=input$p0, dr=input$dp)
        kfn <- modelFun(model=input$degModel, t0=input$startDeg, t1=input$endDeg, 
            tmax=input$Tf, r0=input$k0, dr=input$dk)
        df <- data.frame(t=tseq, p=pfn(tseq), k=kfn(tseq))
        gp <- ggplot(df, aes(t, p))+geom_line(size=1)+xlab("time (minutes)")+
            ylab("production rate (molecules per minute)")+myTheme
        gk <- ggplot(df, aes(t, k))+geom_line(size=1)+xlab("time (minutes)")+
            ylab(expression(paste("degradation rate (minutes"^{-1}, ")")))+myTheme
        pushViewport(viewport(layout=grid.layout(1, 2)))
        print(gp, vp=viewport(layout.pos.row=1, layout.pos.col=1))
        print(gk, vp=viewport(layout.pos.row=1, layout.pos.col=2))
    })
    
    output$solutionPlot <- renderPlot({
        tseq <- seq(0, input$Tf, length=n)
        pfn <- modelFun(model=input$prodModel, t0=input$startProd, t1=input$endProd, 
            tmax=input$Tf, r0=input$p0, dr=input$dp)
        kfn <- modelFun(model=input$degModel, t0=input$startDeg, t1=input$endDeg, 
            tmax=input$Tf, r0=input$k0, dr=input$dk)
        df <- data.frame(t=tseq, p=pfn(tseq), k=kfn(tseq))
        m1 <- log(2)/input$t1; m2 <- log(2)/input$t2
        x01 <- c(X1=x0ss(input$p0, m1, input$k0), 
                 X2=x1ss(input$p0, m1, input$k0))
        x02 <- c(X1=x0ss(input$p0, m2, input$k0), 
                 X2=x1ss(input$p0, m2, input$k0))
        sol <- solveModel(x01=x01, x02=x02, tseq=tseq, m1=m1, m2=m2,
            pfn=pfn, kfn=kfn)
        solMelt <- melt(sol, id="time")
        solMelt <- filter(solMelt, variable %in% c("FP1", "FP2"))
        g1 <- ggplot(solMelt, aes(time, value, color=variable))+geom_line(size=1)+
            scale_color_manual(values=c("darkgreen", "red"))+myTheme+theme(legend.position="none")+
            xlab("time (minutes)")+ylab("FP1, FP2 fluorescence intensity (a.u.)")
        g2 <- ggplot(sol, aes(time, Ratio))+geom_line(color="blue", size=1)+
            scale_color_manual(values=c("darkgreen", "red"))+myTheme+
            xlab("time (minutes)")+ylab("FP2/FP1 intensity ratio")
        pushViewport(viewport(layout=grid.layout(1, 2)))
        print(g1, vp=viewport(layout.pos.row=1, layout.pos.col=1))
        print(g2, vp=viewport(layout.pos.row=1, layout.pos.col=2))
    })
    
    
})
