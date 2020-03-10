#! Noch nicht fertig
step <- function(state,L,Nd,Nu,t,U) {
	N <- Nd + Nu
	#choose electron
	i <- which(state == 1)[sample(1:N, 1)]
	#get neighbors
	neighbors <- c(0,0)
	if(i == 1 | i == L + 1){
		neighbors[1] <- i - 1 + L
	}else {
		neighbors[1] <- i - 1
	}
	if(i == L | i == 2 * L){
		neighbors[2] <- i + 1 - L
	}else{
		neighbors[2] <- i + 1 
	}
	#choose direction
	f <- neighbors[sample(1:2, 1)]
	#check if possible
	if(state[f] == 1) return(FALSE)
	cat(sprintf("%i -> %i", i, f))
	
	#TODO#################################################################################
	#hopping prefactor
	x <- 1.
	if((i %% L == 1 & f %% L == 0) | (f %% L == 1 &  i %% L == 0)){
		if(i > L){
			x <- (-1)^(Nu - 1)
		}else{
			x <- (-1)^(Nd - 1)
		}
	}
	#calculate energy difference
	E <- -U * state[(i - 1 + L) %% (2 * L ) + 1] + U * state[(f - 1 + L) %% (2 * L ) + 1]	#t*???
	p <- exp(-2*E)
	#TODOEND##############################################################################

	if(p >= 1 | runif(1,1) < p){
		state[i] = 0
		state[f] = 1
		return(state)
	}else{
		return(FALSE)
	}
}