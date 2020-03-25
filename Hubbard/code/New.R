M_1 <- function(N_t = 4, phi_t = rep(0, N_t)){
    M <- diag(x=1, nrow=N_t)
    expPhi <- -exp(phi_t)
    #print(expPhi)
    M[row(M) == (col(M) + 1)] <- expPhi[1:(N_t-1)]
    M[1,N_t] <- - expPhi[N_t]
    return(M)
}

Z_1 <- function(N = 50000, N_t = 24, U = 2, beta = 2){
    sum <- 0
    for(i in 1:N){
        phi <- samplePhi(N_t=N_t, U=U, beta=beta)
        sum <- sum + det(M_1(N_t=N_t, phi_t=phi) %*% M_1(N_t=N_t, phi_t=-phi))
    }
    Z <- sum/N
    return(Z)
}

Z_1e <- function(U = 2, beta = 2){
    Z <- 2*(1+exp(beta*U/2))
    return(Z)
}

plot_Z_1_U <- function(beta=2,N=50000){
	U <- list(0,0.5,1,1.5,2)
	data <- c()
    for(V in U){
        data <- c(data,Z_1(U=V, beta=beta))
    }
	pdf(file="plot.pdf")
    plot(x=U,y=data,ylab="Z",main="1-Site Z")
    U <- seq(0,2,0.1)
    lines(U, Z_1e(U=U, beta=beta))
    dev.off()
	return(invisible(matrix(data=c(U,data),ncol=2)))
}

C_1 <- function(N = 50000, N_t = 24, U = 2, beta = 2,tau=0){
    sumC <- 0
    sumZ <- 0
    for(i in 1:N){
        phi <- samplePhi(N_t=N_t, U=U, beta=beta)
        M <- M_1(N_t=N_t, phi_t=phi)
        det <- det(M %*% M_1(N_t=N_t, phi_t=-phi))
        sumC <- sumC + det * solve(M)[,1]
        sumZ <- sumZ + det
    }
    C <- sumC/sumZ
    return(C)
}

C_1e <- function(U=2, beta = 2 ,tau=0){
    C <- cosh(U/2*(tau-beta/2))/(2*cosh(U*beta/4))
    return(C)
}

plot_C_1_tau <- function(beta=2, U=2, N=50000){
	N_t <- 24
    tau <- seq(0,beta-beta/N_t, length.out= N_t)
    data <- C_1(tau=tau, U=U, beta=beta, N_t=N_t)
    pdf(file="plot.pdf")
    plot(x=tau, y=data, ylab="C", main="1-Site C")
    tau2 <- seq(0,2,0.1)
    lines(tau,C_1e(U=U, beta=beta, tau=tau))
    dev.off()
	return(invisible(matrix(data=c(tau,data),ncol=2)))
}

samplePhi <- function(N_t = 64, U = 5, beta = 3){
    U <- U*beta/N_t
    phi <- rnorm(n=N_t,mean=0,sd=sqrt(U))
    return(phi)
}

Figure1a <- function(N_t = 32){
    results <- eigen(M_1(N_t=N_t, phi_t = rep(0,N_t)))$values
    plot(results)
}

Figure1b <- function(sets =20, N_t = 64,U = 5, beta=3){
    results <- c()
    for(x in 1:20){
        results <- c(results,eigen(M_1(N_t=N_t, phi_t = rnorm(n=N_t,mean=0,sd=sqrt(U*beta/N_t))))$values)
    }
    plot(results)
}
