expE <- function(beta=1, U=3) {
	En <- c(0,U,(U+sqrt(16+U*U))/2,(U-sqrt(16+U*U))/2)
	expE <- exp(-beta*En) %*% En / sum(exp(-beta*En))
	return(expE[1,1])
}

plot_beta <- function(U=3){
	beta <- list(0,0.1,1,10)
	plot(x=beta,y=lapply(beta, expE, U=U))
}

plot_U <- function(beta=1){
	U <- list(0,0.5,1,2,4,8,16)
	plot(x=U,y=lapply(U, expE, beta=beta))
}