expE <- function(beta=1, U=3) {
	En <- c(0,U,(U+sqrt(16+U*U))/2,(U-sqrt(16+U*U))/2)
	expE <- exp(-beta*En) %*% En / sum(exp(-beta*En))
	return(expE[1,1])
}

expProp <- function(beta=1, U=3,i=1,j=2){
	t <- 1
	En <- c(0,U,(U+sqrt(16+U*U))/2,(U-sqrt(16+U*U))/2)
	eig <- list(c(0,-1,1,0),c(1,0,0,-1),c(1,En[4]/(2*t),En[4]/(2*t),1),c(1,En[3]/(2*t),En[3]/(2*t),1))
	expE <- 0
	for(n in 1:4){
		state <- eig[[n]]
		expE <- expE + exp(-beta*En[n])* (state %*% Prop(i=i,j=j) %*% state) / (state %*% state)
	}
	expE <- expE / sum(exp(-beta*En))
	return(expE)
}

Prop <- function(i=1,j=2){
	Basis <- list(c(0,1,0,1),c(0,1,1,0),c(1,0,0,1),c(1,0,1,0))
	m <- matrix(data=rep(0,16),nrow=4)
	for(state1 in 1:4){
		initial <- Basis[[state1]]
		if(initial[i]==1&(initial[j]==0|i==j)){
			x <- initial
			x[i] <- 0
			x[j] <- 1
			for(state2 in 1:4){
				final <- Basis[[state2]]
				if(all(x==final)) m[state1,state2] <- 1
			}
		}
	}
	return(m)
}

plot_beta <- function(U=3){
	beta <- list(0,0.1,0.5,1,10)
	plot(x=beta,y=lapply(beta, expE, U=U))
}

plot_U <- function(beta=1){
	U <- list(0,0.5,1,2,4,8,16)
	plot(x=U,y=lapply(U, expE, beta=beta))
}