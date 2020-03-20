expect <- function(Psi, Basis, Operator) {
	sum_n <- 0
	sum_d <- 0
	for(state in Basis){
		xPsi <- Psi[which(Psi$state==state2int(state)), ]$factor
		sum_d <- sum_d + xPsi^2
		sum_n <- sum_n + xPsi^2*OP_loc(state, Operator, Psi)
	}
	return(sum_n/sum_d)
}

OP_loc <- function(x, O, Psi){
	xOPsi <- 0
	x <- data.frame(state = c(state2int(x)), factor = c(1.0))
	y <- O(x)		#dataframe of reached basis states 
	xPsi <- Psi[which(Psi$state==state2int(x)), ]$factor
}

state2int <- function(state){
	int <- 0
	pos <- 0
	for(x in state){
		int <- int + x*2^pos
		pos <- pos + 1
	}
	return(int)
}

basis2vec <- function(basis){
	vec <- c()
	count <- 1
	for(state in basis){
		vec[count] <- state2int(state)
		count <- count + 1
	}
	return(vec)
}

TestState <- c(0,1,0,1)
TestBasis <- list(c(0,1,0,1), c(1,0,0,1), c(0,1,1,0), c(1,0,1,0))
TestPsi <- data.frame(state = basis2vec(TestBasis), factor = c(1,2,3,4))