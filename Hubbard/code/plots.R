source("New.R", encoding = "UTF-8")

plot_Z1_Nt <-function(N=50000, Nt=c(4,8,16,32,64), U=2, beta=2){
	x <- Nt
    y <- c()
	dy <- c()
    for(n in Nt){
        Z <- Z(N=N, L=1, N_t=n, U=U, beta=beta)
        y <- c(y,Z$mean)
        dy <- c(dy, Z$sd)
    }
    out <- data.frame(Nt=x,Z=y,dZ=dy,N,beta,U)
    write.csv(out, "plot_data.csv")
    Ze <- Z_1e(U=U, beta=beta)
    pdf(file="plot.pdf")
    plotwitherror(x=x,y=y,dy=dy,xlab="N_t",ylab="Z",main="1-Site U=2 beta=2 N=50000",xlim=c(0,tail(Nt,n=1)+1),ylim=c(Ze-2, Ze+2),col="blue")
    #lines(Nt[1], Ze)
    abline(h=Ze,col="red")
    legend("topright",legend=c("MonteCarlo results with errors","exact result"),bty="n",col=c("blue","red"),pt.bg="blue",pch=c(21,19))
    dev.off()
	return(invisible(out))
}

plot_Z1_N <-function(N=c(1000,5000,10000,20000,50000,100000), Nt=32, U=2, beta=2){
	x <- N
    y <- c()
	dy <- c()
    for(n in N){
        Z <- Z(N=n, L=1, N_t=Nt, U=U, beta=beta)
        y <- c(y,Z$mean)
        dy <- c(dy, Z$sd)
    }
    out <- data.frame(N=x,Z=y,dZ=dy,Nt,beta,U)
    write.csv(out, "plot_Z1N_data.csv")
    Ze <- Z_1e(U=U, beta=beta)
    pdf(file="plot_Z1N.pdf")
    plotwitherror(x=x,y=y,dy=dy,xlab="N",ylab="Z",main="1-Site U=2 beta=2 Nt=32",xlim=c(0,tail(N,n=1)+1000),ylim=c(Ze-4, Ze+4),col="blue")
    abline(h=Ze,col="red")
    legend("topright",legend=c("MonteCarlo results with errors","exact result"),bty="n",col=c("blue","red"),pt.bg="blue",pch=c(21,19))
    dev.off()
	return(invisible(out))
}

plot_Z1_U <- function(N=50000, Nt=32, U=c(0,1,2,4,8), beta=2){
	x <- U
    y <- c()
	dy <- c()
    for(u in U){
        Z <- Z(N=N, L=1, N_t=Nt, U=u, beta=beta)
        y <- c(y,Z$mean)
        dy <- c(dy, Z$sd)
    }
    out <- data.frame(U=x,Z=y,dZ=dy,Nt,beta,N)
    write.csv(out, "plot_Z1U_data.csv")
    pdf(file="plot_Z1U.pdf")
    plotwitherror(x=x,y=y,dy=dy,xlab="U",ylab="Z",main="1-Site N=50000 beta=2 Nt=32",xlim=c(0,tail(U,n=1)+1),ylim=c(0,tail(y+dy,n=1)+100),col="blue")
    U <- seq(0,8,0.1)
    lines(U, Z_1e(U=U, beta=beta),col="red")
    legend("topleft",legend=c("MonteCarlo results with errors","exact result"),bty="n",col=c("blue","red"),pt.bg="blue",pch=c(21,19))
    dev.off()
	return(invisible(out))
}

plot_Z1_b <- function(N=50000, Nt=32, U=2, beta=c(0,1,2,4,6)){
	x <- beta
    y <- c()
	dy <- c()
    for(b in beta){
        Z <- Z(N=N, L=1, N_t=Nt, U=U, beta=b)
        y <- c(y,Z$mean)
        dy <- c(dy, Z$sd)
    }
    out <- data.frame(beta=x,Z=y,dZ=dy,Nt,U,N)
    write.csv(out, "plot_Z1b_data.csv")
    pdf(file="plot_Z1b.pdf")
    plotwitherror(x=x,y=y,dy=dy,xlab="beta",ylab="Z",main="1-Site N=50000 U=2 Nt=32",xlim=c(0,tail(beta,n=1)+1),col="blue")
    beta <- seq(0,8,0.1)
    lines(beta, Z_1e(U=U, beta=beta),col="red")
    legend("topleft",legend=c("MonteCarlo results with errors","exact result"),bty="n",col=c("blue","red"),pt.bg="blue",pch=c(21,19))
    dev.off()
	return(invisible(out))
}

plot_C2_t <- function(N=10000,Nt=32,U=2,beta=2){
    x <- seq(0,beta*(Nt-1)/Nt,length.out=Nt)
    y <- c()
	dy <- c()
    C <- C(N=N, N_t=Nt,L=2, U=U, beta=b)
    y <- C$mean
    dy <- C$sd
    out <- data.frame(tau=x,C=y,dC=dy,Nt,U,N)
    write.csv(out, "plot_C2t_data.csv")
    pdf(file="plot_C2t.pdf")
    plotwitherror(x=x,y=y[1:Nt],dy=dy[1:Nt],xlab="tau",ylab="Correlator",main="2-Site U=2 beta=2",xlim=c(0,beta),ylim=c(0.9*min(y)-0.05,1.1*max(y)))
    pointswitherror(x=x,y=y[(Nt+1):(2*Nt)],dy=dy[(Nt+1):(2*Nt)])
    #beta <- seq(0,8,0.1)
    #lines(beta, Z_1e(U=U, beta=beta),col="red")
    legend("topright",legend=c("C_11","C_12"),bty="n",col=c("blue","green"),pt.bg="blue",pch=1)
    dev.off()
	return(invisible(out))
}

plotwitherror<-function(x,y,dy,col="blue",...) {
    plot(x=x,y=y,col=col, ...)
    arrows(x0=x,y0=y-dy,x1=x,y1=y+dy,length=0.01,angle=90,code=3,col=col)
}

pointswitherror<-function(x,y,dy,col="green",...) {
    points(x=x,y=y,col=col, ...)
    arrows(x0=x,y0=y-dy,x1=x,y1=y+dy,length=0.01,angle=90,code=3,col=col)
}

#out <- data.frame(x,y,dy)
#write.csv(out, "plot_data.csv")