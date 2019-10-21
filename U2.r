plotwitherror<-function(x,y,dy,col="black",...) {
    plot(x=x,y=y,col=col, ...)
    arrows(x0=x,y0=y-dy,x1=x,y1=y+dy,length=0.01,angle=90,code=3,col=col)
}

Nobs = 1000
mu <- rep(0,25)
var <- rep(0,25)
for (n in 1:25) {
    z <- rep(0,Nobs)
    for (obs in 1:Nobs) {
        for (i in 1:n) {
            z[obs] <- z[obs] + runif(1)
        }
        z[obs] <- z[obs] - n/2
        z[obs] <- z[obs] * (12/n)**0.5
    }
    mu[n] <- mean(z)#sum(z)/Nobs
    var[n] <- var(z)
}
plotwitherror(1:25,mu,var,col="blue",xlim=c(1,25),ylim=c(-1.5,1.5))