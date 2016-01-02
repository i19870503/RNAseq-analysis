gap.barplot <- function(data, btm, top, ratio=1, sd=NULL, 
                        brk.size=1.5, brk.srt=145, brk.col='black', brk.bg='white', ...){
    x <- rownames(data)
    y <- data
	#y <- t(data[,- 1])
    if(missing(sd)) sdx=0 else sdx <- t(sd)
    sdu <- y + sdx
    sdd <- y - sdx
    ylim <- c(0, max(sdu) * 1.05)
    ## Set up virtual y limits
    halflen <- btm - ylim[1]
    xlen <- halflen * 0.1
    v_tps1 <- btm + xlen        # virtual top positions
    v_tps2 <- v_tps1 + halflen * ratio
    v_ylim <- c(ylim[1], v_tps2)
    r_tps1 <- top               # real top positions
    r_tps2 <- ylim[2]
    ## Rescale data
    lmx <- summary(lm(c(v_tps1, v_tps2)~c(r_tps1, r_tps2)))
    lmx <- lmx$coefficients
    sel1 <- y > top
    sel2 <- y >=btm & y <=top
    y[sel1] <- y[sel1] * lmx[2] + lmx[1]
    y[sel2] <- btm + xlen/2

    sel1 <- sdd > top
    sel2 <- sdd >=btm & sdd <=top
    sdd[sel1] <- sdd[sel1] * lmx[2] + lmx[1]
    sdd[sel2] <- btm + xlen/2

    sel1 <- sdu > top
    sel2 <- sdu >=btm & sdu <=top
    sdu[sel1] <- sdu[sel1] * lmx[2] + lmx[1]
    sdu[sel2] <- btm + xlen/2
    ## bar plot
    xx <- barplot2(y, beside=1, ylim=v_ylim, axes = FALSE, col=c('blue','red') ,xaxt='n',ylab='Relative Expression',main='RT-qPCR results of piPS vs pEF',...)
    xps <- 8
    n <- nrow(y)
    #if( n > 1) xps <- colMeans(xx)
    #axis(1, at=8, label=x, ...)
    ## error bars
    if(!missing(sd)){
        bw <- 0.2
        for(i in 1:n){
            segments(xx[i,], sdu[i,], xx[i,], sdd[i,], lwd=2)
            segments(xx[i,]-bw, sdu[i,], xx[i,]+bw, sdu[i,], lend=2)
            segments(xx[i,]-bw, sdd[i,], xx[i,]+bw, sdd[i,], lend=2)
        }
    }
    ## Real ticks and labels    
    brks1 <- pretty(seq(0, btm, length=10), n=4)
    brks1 <- brks1[brks1 >= 0 & brks1 < btm]
    brks2 <- pretty(seq(top, r_tps2, length=10), n=4)
    brks2 <- brks2[brks2 > top & brks2 <= r_tps2]
    labx <- c(brks1, brks2)
    ## Virtual ticks
    brks <- c(brks1, brks2 * lmx[2] + lmx[1])
    axis(2, at=brks, labels=labx, ...)
    box()
    ## break marks
    pos <- par("usr")[1:2]
    pox <- (pos[2] - pos[1])/100
    rect(pos[1] - pox, btm, pos[2] + pox, v_tps1, col=brk.bg, xpd=TRUE, border=brk.bg)
    xp <- c(pos[1], pos[1], pos[2], pos[2])
    yp <- c(btm, v_tps1, btm, v_tps1)
    text(xp, yp, expression(NULL - NULL), cex=brk.size, xpd=T, adj=c(0.5, 1), srt=brk.srt)
    invisible(xx)
}
