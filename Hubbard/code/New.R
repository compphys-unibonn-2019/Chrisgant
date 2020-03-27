H2 <- function(U=2){
    B <- matrix(data=c(0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1),nrow=4)
    H <- diag(x=0, nrow=16)
    for(i in 1:16){
        si <- B[,i]
        nui <- sum(si[c(1,2)])
        ndi <- sum(si[c(3,4)])
        for(j in 1:16){
            sj <- B[,j]
            if(i==j){
                H[i,j] <- -U/2*((si[1]-si[3])^2+(si[2]-si[4])^2)
            }else{
                nuj <- sum(sj[c(1,2)])
                ndj <- sum(sj[c(3,4)])
                if(nui == 1 & nuj == 1 & ndi == ndj & si[3]==sj[3]) H[i,j] <- -1
                if(ndi == 1 & ndj == 1 & nui == nuj & si[1]==sj[1]) H[i,j] <- -1
            }
            
        }
    }
    return(H)
}

M_1 <- function(N_t = 4, phi_t = rep(0, N_t)){
    M <- diag(x=1, nrow=N_t)
    expPhi <- -exp(phi_t)
    #print(expPhi)
    M[row(M) == (col(M) + 1)] <- expPhi[1:(N_t-1)]
    M[1,N_t] <- - expPhi[N_t]
    return(M)
}

M <- function(N_t = 4, N = 2, phi_t = rep(0, N_t*N), t=1,beta=2){  #! not ready for continuous boundary
    M <- diag(x=1, nrow=N_t*N)
    t <- t*beta/N_t
    expPhi <- exp(-phi_t)
    M[row(M) == (col(M) + 1) & col(M) %% N_t != 0]   <- -expPhi[1:(N*N_t-N)]
    M[row(M) == (col(M) - N_t + 1) & col(M) %% N_t == 0]   <- +expPhi[(N*N_t-N+1):(N_t*N)]
    M[row(M) -1 == (col(M)) %% (2*N_t) & col(M) %% N_t == 0]   <- -t
    M[row(M) -1 == (col(M) - N_t) %% (2*N_t) & col(M) %% N_t != 0]   <- +t
    return(M)
}

Z_1 <- function(N = 50000, N_t = 24, U = 2, beta = 2){
    sum <- 0
    MM <- rep(0,N)
    n <- 0
    for(i in 1:N){
        phi <- samplePhi(N_t=N_t, U=U, beta=beta)
        x <- det(M_1(N_t=N_t, phi_t=phi) %*% M_1(N_t=N_t, phi_t=-phi))
        MM[i] <- x
        sum <- sum + x
        if(x < 0) n <- n +1
    }
    #print(n)
    #print(max(MM))
    #hist(MM[which(MM < 100)], breaks=100)
    Z <- sum/N
    return(Z)
}

Z <- function(N = 50000, N_t = 4, L = 2, U = 2, beta = 2){
    sum <- 0
    MM <- rep(0,N)
    n <- 0
    for(i in 1:N){
        phi <- samplePhi(N_t=N_t, N=L, U=U, beta=beta)
        x <- det(M(N_t=N_t, N=L, phi_t=phi,beta=beta) %*% M(N_t=N_t, N=L, phi_t=-phi,beta=beta))
        MM[i] <- x
        sum <- sum + x
        if(x < 0) n <- n +1
    }
    #print(n)
    hist(MM[which(MM < 3000)], breaks=80)
    Z <- sum/N
    return(Z)
}

Z_1e <- function(U = 2, beta = 2){
    Z <- 2*(1+exp(beta*U/2))
    return(Z)
}

Z_2e <- function(U = 2, beta = 2, t=1){
    Z <- 3+3*exp(beta*U)+2*exp(beta*U/2)*(4*cosh(beta*t)+cosh(beta*sqrt(U^2/4+4*t^2)))
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

plot_C_1_Nt <- function(beta=2, U=2, N=50000,N_t=24){
	#N_t <- 48
    tau <- seq(0,beta-beta/N_t, length.out= N_t)
    data <- C_1(tau=tau, U=U, beta=beta, N_t=N_t)
    pdf(file="plot.pdf")
    plot(x=tau, y=data, ylab="C", main="1-Site C")
    tau2 <- seq(0,2,0.1)
    lines(tau,C_1e(U=U, beta=beta, tau=tau))
    dev.off()
	return(invisible(matrix(data=c(tau,data),ncol=2)))
}

samplePhi <- function(N_t = 64, U = 5, beta = 3,N=1){
    U <- U*beta/N_t
    phi <- rnorm(n=N_t*N,mean=0,sd=sqrt(U))
    #phi <- matrix(data=phi,ncol=N)
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
