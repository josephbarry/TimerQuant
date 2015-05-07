plotPrimordiumProfile <- function(x, add=FALSE, ylab="ratio", lwd=2, 
    cex.lab=1.5, cex.axis, xlim=NA, ylim=NA, main="", col="black", lty=1, 
    alpha=0.3) {
    xPts <- as.numeric(colnames(x))
    yMedian <- apply(x, 2, median, na.rm=TRUE)
    if (sum(is.na(yMedian)) == length(yMedian)) return(0)
    yMad <- apply(x, 2, mad, na.rm=TRUE)
    xlab <- "distance from leading edge (microns)"
    yPoly <- c(yMedian+yMad, rev(yMedian-yMad))
    yPoly[is.na(yPoly)] <- 0
    xPoly <- -c(xPts, rev(xPts))
    if (is.na(ylim[1])) ylim <- c(0, max(yPoly, na.rm=TRUE))
    if (is.na(xlim[1])) xlim <- c(min(-xPts, na.rm=TRUE), 0)
    if(!add) {
        par(mar=c(5, 5, 4, 1)+0.1)
        plot(-xPts, yMedian, xlab=xlab, ylab=ylab, type="l", lwd=lwd, 
            axes="F", cex.lab=cex.lab, xlim=xlim, ylim=ylim, main=main, col=col,
            lty=lty)
        axisSeq <- seq(-200, 0, by=25)
        axis(1, at=axisSeq, labels=-axisSeq, cex.axis=cex.axis)
        axis(2, cex.axis=cex.axis)
    }
    else {
        points(-xPts, yMedian, type="l", lwd=lwd, col=col, lty=lty)
    }
    polygon(xPoly, yPoly, col=rgb(red=0, green=0, blue=0.5, alpha=alpha), 
        border=NA)
}
