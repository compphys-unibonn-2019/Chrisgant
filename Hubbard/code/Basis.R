#Returns basis for different conditions: Lattice site set, Lattice site and Number of e set, Lattice site and Number of spin up and down set
install.packages('gtools')
library(gtools)

L_set <- function(L) {
    return(permutations(n=2,r=2*L,v=0:1,repeats.allowed=TRUE))
}

N_set <- function(L, N) {
    x <- L_set(L)
    s <- rowSums(x)
    return(x[which(s == N),])
}

Nd_set <- function(L, Nd, Nu){
    x <- N_set(L, Nu+Nd)
    if(is.null(nrow(x))){
        s <- sum(x[1:L])
        if(s == Nd){
            return(x)
        }else {
            return(NULL)
        }
    }
    s <- rowSums(x[,1:L])
    return(x[which(s == Nd),])
}

#############Test area#############################
#Nd_set(5,2,1)