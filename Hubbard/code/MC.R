source("Basis.R")
source("Hamiltonian.R")

N <- 10
beta <- 0.1

HubbardInit <- function(L, Nd = L/2, Nu = L/2){
	S <- array(rep(0, 2*L), dim=2*L)
	S[sample(1:L,Nd)] <- 1
	S[sample((L+1):(2*L),Nu)] <- 1
	return(S)
}

deltaH <- function(h, S, t=1, U=3, continuous=TRUE){	#! unsure				H|Psi> = E|Psi'> => E=|Psi'|
	L <- length(S)/2
	dH_pot <- -U * S[(h[1] - 1 + L) %% (2 * L ) + 1] + U * S[(h[2] - 1 + L) %% (2 * L ) + 1]
	#neighbors_i <- getNeighbors(L=L, i=h[1], continuous=continuous)
	#neighbors_f <- getNeighbors(L=L, i=h[2], continuous=continuous)
	dH_kin <- -t# * (2 * sum(S[neighbors_i[!is.na(neighbors_i)]]) - 2 * (sum(S[neighbors_f[!is.na(neighbors_f)]]) - 1))
	dH <- dH_pot + dH_kin
	return(dH)
	#? Energy expectation value?
}

H <- function(S, t=1, U=3, continuous=TRUE){
	L <- length(S)/2
	H_pot <- U * (S[1:L] %*% S[(L+1):(2*L)])[1]
	n <- 0
	for(i in which(S==0)){
		n <- n + sum(S[getNeighbors(L=L,i=i,continuous=continuous)])
	}
	H_kin <- -t * n
	H <- H_pot + H_kin
	return(H)
}

getNeighbors <- function(L, i, continuous=TRUE){
	neighbors <- c()
	if(i == 1 | i == L + 1){
		if(continuous) neighbors[1] <- i - 1 + L
	}else {
		neighbors[1] <- i - 1
	}
	if(i == L | i == 2 * L){
		if(continuous) neighbors[2] <- i + 1 - L
	}else{
		neighbors[2] <- i + 1
	}
	return(neighbors)		# neighbors[1] site links, neighbors[2] site rechts
}

HubbardSweep <- function(S, t=1, U=3, beta=1,continuous=TRUE){
	L <- dim(S)/2
	AR <- 0
	rand <- runif(n=2*L)
	for(x in sample(1:(2*L), 2*L)){					#? maybe random order?  x in sample(1:(2*L), 2*L)
		y <- getNeighbors(L,x,continuous=continuous)[2]
		if(!is.na(y) & (S[x]+S[y] == 1)){
			i <- x
			f <- y
			if(S[i] == 0){
				i <- y
				f <- x
			}
			dH <- deltaH(h=c(i,f), S=S, t=t, U=U)
			Accept <- (rand[x] < exp(-beta*dH))
			if(Accept){
				S[i] = 0
				S[f] = 1
				AR <- AR+1
			}
		}
	}
	return(list(S=S, AcceptanceRate=AR/(2*L)))
}

HubbardPlot <- function(S){
	cc=c("white","blue")
	L <- dim(S)/2
	df <- data.frame(down = S[1:L], up = S[(1+L):(2*L)])
	#library(ggplot2)
	#df <- data.frame(x = 1:L, y = S[1:L]+S[(1+L):(2*L)])
	#ggplot(df, aes(x=x, y=0, fill = y)) + geom_tile()
	heatmap(x=t(as.matrix(df)), Rowv = NA, Colv = NA, scale="none", col=cc, xlab="site", ylab="down   up   ", labRow=NA, labCol=NA, margins=c(0.5,0.5), main="Configuration")
}

Metropolis <- function(N=100, L=6, Nd = L/2, Nu = L/2, beta=1, t=1, U=3,continuous=TRUE){
	Basis <- GenerateBasis(N=N, L=L, Nd=Nd, Nu=Nu, beta=beta, t=t, U=U, continuous=continuous)
	for(i in 1:nrow(Basis)){
		E[i] <- H(S)
	}
	print(sum(E)/N)
	return(invisible(S))
}

GenerateBasis <- function(N=100, restarts=1, L=6, Nd = L/2, Nu = L/2, beta=1,max=2000, t=1, U=3,continuous=TRUE){
	S <- HubbardInit(L=L, Nd=Nd, Nu=Nu)
	Basis <- matrix(data=S,nrow=1,ncol=2*L)
	size <- 0
	for(r in 1:restarts){
		if(r > 1)S <- HubbardInit(L=L, Nd=Nd, Nu=Nu)
		#if(any(apply(Basis, 1, function(x, want) isTRUE(all.equal(x, want)), c(S))) == FALSE){
			Basis <- rbind(Basis,t(S))
			size <- size + 1
		#}
		for(i in 2:N){
			#print(S)
			S <- HubbardSweep(S=S, beta=beta, t=t, U=U, continuous=continuous)$S
			#if(any(apply(Basis, 1, function(x, want) isTRUE(all.equal(x, want)), c(S))) == FALSE){
				Basis <- rbind(Basis,t(S))
				size <- size + 1
			#}
			if(size == max) break
		}
	}
	
	#print(sum(E)/N)
	return(Basis)
}

H_Sub <- function(Basis,t=1,U=3) {
	L <- ncol(Basis)
	n <- nrow(Basis)
	if(is.null(n)){
		if(length(Basis) > 0){
			n <- 1
			b <- matrix(data = Basis, nrow = 1)
		}else {
			return(NULL)
		}
	}
	if(n == 0) return(NULL)
	H <- matrix(data = 0, nrow = n, ncol = n)
	for (i in 1:n) {
		for (j in 1:i) {
			x <- H1(Basis[i,],Basis[j,],t,U,TRUE)
			H[i,j] = x
			H[j,i] = x
		}
	}
	return(H)
}