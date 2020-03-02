source("Basis.R")
source("Hamiltonian.R")

H_L <- function(L,t,U) {
    b <- L_set(L)
    n <- nrow(b)
    #print(n)
    H <- matrix(data = 0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:i) {
            x <- H1(b[i,],b[j,],t,U,TRUE)
            H[i,j] = x
            H[j,i] = x
        }
    }
    return(H)
}
#L_set(1)
#H_L(1,1,3)

H_N <- function(L,N,t,U) {
    b <- N_set(L,N)
    n <- nrow(b)
    #print(n)
    H <- matrix(data = 0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:i) {
            x <- H1(b[i,],b[j,],t,U,TRUE)
            H[i,j] = x
            H[j,i] = x
        }
    }
    return(H)
}
#N_set(2,2)
#H_N(2,2,1,3)

H_Nd <- function(L,Nd,Nu,t,U) {
    b <- Nd_set(L,Nd,Nu)
    n <- nrow(b)
    if(is.null(n)){
        if(length(b) > 0){
            n <- 1
            b <- matrix(data = b, nrow = 1)
        }else {
            return(NULL)
        }
    }
    if(n == 0) return(NULL)
    H <- matrix(data = 0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:i) {
            x <- H1(b[i,],b[j,],t,U,TRUE)
            H[i,j] = x
            H[j,i] = x
        }
    }
    return(H)
}

H_L_split <- function(L,t,U){
    for(N in 0:(2*L)){
        for(Nd in 0:L){
            for(Nu in 0:L){
                if(Nd+Nu==N){
                    cat(sprintf("Nd = %i  Nu = %i: \n", Nd, Nu))
                    print(H_Nd(L,Nd,Nu,t,U))
                }
            }
        }
    }
}
#TODO return vector of matrices instead of print

H_N_split <- function(L,N,t,U){    
    for(Nd in 0:N){
        for(Nu in 0:N){
            if(Nd+Nu==N){
                cat(sprintf("Nd = %i  Nu = %i: \n", Nd, Nu))
                print(H_Nd(L,Nd,Nu,t,U))
            }
        }
    }
}
#TODO return vector of matrices instead of print

#############Test area#############################
#Nd_set(2,0,0)
#H_Nd(2,1,1,1,3)
#H_N_split(2,2,1,3)


eigen(H_L(1,1,3))