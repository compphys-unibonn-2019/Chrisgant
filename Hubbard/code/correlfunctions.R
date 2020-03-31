plot_C11_analytic(beta=2,U=2,t=1,tau){
    w <- sqrt(4*t^2+U^2/4)
    alpha <-0.37
    beta <- 0.60
    Z <- 3*(1+exp(beta*U))+2*exp(beta*U/2)*(4*cosh(beta*t)+cosh(beta*w)
    C11 <- (1/Z)*(1+exp(beta*U)*(exp(tau*(U/2-w))+(1/2)*exp(tau*(U/2-t)))+(1/2)*exp(beta*(U/2-t))*(exp(tau*(w-t))+2+exp(tau*(-U/2-t)))+(1/2)*exp(beta*(U/2+t))*(exp(tau*(w-t))+exp(tau*(-U/2+t)))+exp(beta*(U/2+w))*exp(tau*(t-w))*(alpha^2+beta^2)+exp(beta*(U/2-w))*exp(tau*(-w-t))*(alpha^2+beta^2))
    pdf(file="plot.pdf")

    plot(x=tau,y=C11,xlab="tau",ylab="C_11",main="2-Site U=2 beta=2 analytic",xlim=c(0,2),ylim=c(0, 2),col="blue")
    legend("topright",legend=c("MonteCarlo results with errors","exact result"),bty="n",col=c("blue","red"),pt.bg="blue",pch=c(21,19))
    dev.off()
}
