deltaH <- function(h, S, t=1, U=3, continuous=TRUE){
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
	#n <- 0
	#for(i in which(S==0)){
	#	n <- n + sum(S[getNeighbors(L=L,i=i,continuous=continuous)])
	#}
	H_kin <- -t #* n
	H <- H_pot + H_kin
	return(H)
}

HubbardInit <- function(L, Nd = L/2, Nu = L/2){
	S <- array(rep(0, 2*L), dim=2*L)
	S[sample(1:L,Nd)] <- 1
	S[sample((L+1):(2*L),Nu)] <- 1
	return(S)
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

Metropolis <- function(N=10000, restarts=1, L=6, Nd = L/2, Nu = L/2, beta=1, t=1, U=3,continuous=TRUE,observable=H,...){
	S <- HubbardInit(L=L, Nd=Nd, Nu=Nu)
	E <- numeric()
	for(r in 0:(restarts-1)){
			for(i in 1:N){
			S <- HubbardSweep(S=S, beta=beta, t=t, U=U, continuous=continuous)$S
			E[r*N+i] <- observable(S)
		}
	}
	return(sum(E)/(N*restarts))
}