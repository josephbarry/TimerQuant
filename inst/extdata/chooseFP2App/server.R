library(TimerQuant)
library(dplyr)
library(ggplot2)
library(gridExtra)

myTheme <- theme_bw()+theme(title=element_text(size=18), text=element_text(size=12), 
    axis.title.y=element_text(vjust=1.5), plot.title=element_text(size=12))

shinyServer(function(input, output, session) {
    output$timerSignal <- renderPlot({
        T2 <- getSpacedSeq(c(input$t2range[1], input$t2range[2]), n=100)
        res <- lapply(seq_along(T2), function(i)  simulatedSignalN(T1=input$t1, 
            T2=T2[i], TA=input$tA, TB=input$tB, sigmaAdd=input$sigma, 
            N=input$nRealizations, p=input$p, E=input$E))
        df <- lapply(seq_along(res), function(i) data.frame(p=input$p, T=T2[i], D=res[[i]]))
        df <- do.call("rbind", df)
        dfs <- df %>% group_by(T) %>% 
            summarise(D.mean=mean(D, na.rm=TRUE), D.sd=sd(D, na.rm=TRUE))
        dfs$D0 <- simulatedSignal(input$t1, dfs$T, input$tA, 
            input$tB, sigmaAdd=0, p=input$p, E=input$E)
        g1 <- ggplot(dfs, aes(T, D.mean, colour="red"))+geom_point()+myTheme+
            xlab("FP2 maturation time (minutes)")+ylab("timer signal S")+
            scale_x_log10(breaks=getBreaks10(input$t2range))+
            geom_errorbar(aes(ymin=D.mean-D.sd, ymax=D.mean+D.sd))+
            geom_line(aes(T, D0), colour="black", data=dfs)+
            theme(legend.position="none")+
            annotation_logticks(sides="b")
        dfs$CV <- dfs$D.sd/dfs$D.mean
        dfsCV <- fitCV(dfs)
        g2 <- ggplot(dfs, aes(T, CV), color="blue")+geom_point()+myTheme+
            scale_x_log10(breaks=getBreaks10(range(dfsCV$T)))+
            scale_y_log10(breaks=getBreaks10(range(dfsCV$CV)))+
            xlab("FP2 maturation time (minutes)")+ylab(expression(sigma[S]~"/"~mu[S]))+
            geom_line(aes(T, CV), data=dfsCV)+
            geom_vline(aes(xintercept=FP2optimumTime), data=dfsCV, 
            linetype="dashed")+theme(legend.position="none")+
            ggtitle(sprintf("Optimal FP2 maturation time is %0.1f minutes", dfsCV$FP2optimumTime[1]))+
            annotation_logticks()
        g <- grid.arrange(g1, g2)
        return(g)
    })
})
