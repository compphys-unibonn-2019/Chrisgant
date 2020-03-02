source("Basis.R")
source("Hamiltonian.R")

H_L <- function(L,t,U) {
    b <- L_set(L)
    n <- nrow(b)
    print(n)
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
    print(n)
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
    print(n)
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
Nd_set(2,2,0)
H_Nd(2,2,0,1,3)