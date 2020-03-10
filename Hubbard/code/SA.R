h <- function(H, phi) {
    return((phi %*% H %*% phi)/(phi %*% phi))
}

SA <- function(H, phi_0 = (rep(1/sqrt(ncol(H)),ncol(H))), T_0 = 100, alpha = 0.9, max = 300, tries = 10, epsilon = 0.0001){
    #initialise
    phi <- phi_0
    T <- T_0
    h_old <- h(H,phi)
    h_min <- h_old
    phi_min <- phi
    #iterate
    for(i in 1:tries){
        for(j in 1:max){
            rand <- rnorm(length(phi), mean = 0, sd = 0.2)*2                    #
            rand2 <- runif(length(phi),1)                                       #~5%
            for(k in 1:length(phi)){
                phi_new <- phi
                phi_new[k] <- phi_new[k] + rand[k]
                #phi_new <- phi_new / sqrt(sum(phi_new^2))
                #h_new <- h(H,phi_new)											#~70%
				h_new <- (phi_new %*% H %*% phi_new)/(phi_new %*% phi_new)
                accept <- FALSE                                                 
                if(h_new < h_old) accept <- TRUE
                else if(rand2[k] < exp((h_old - h_new)/T)) accept <- TRUE
                if(accept){
                    phi[k] <- phi_new[k]
                    #phi <- phi / sqrt(sum(phi_new^2))
                    h_old <- h_new
                    if(h_new < h_min){
                        h_min <- h_new
                        phi_min <- phi / sqrt(sum(phi^2))
                    }
                }
            }
            T <- T*alpha
        }
    }
    #print(phi_min)
    #print(h_min)
    return(phi_min)
}


source("HMatrix.R")
H <- H_Nd(4,2,2,1,3)

phi = (rep(1/sqrt(ncol(H)),ncol(H)))
h_old = 3
h_min = 3
T = 100

phi_min <- SA(H)
system.time(SA(B))
system.time({
    for(i in 1:(3000)){
        rand <- rnorm(36, mean = 0, sd = 0.2)*2
        rand2 <- runif(36,1)
        for(k in 1:36){
            phi_new <- phi
            phi_new[k] <- phi_new[k] + rand[k]
            #phi_new <- phi_new / sqrt(sum(phi_new^2))
            h_new <- h(H,phi_new)
            accept <- FALSE
            if(h_new < h_old) accept <- TRUE
            else if(rand2[k] < exp((h_old - h_new)/T)) accept <- TRUE
            if(accept){
                phi[k] <- phi_new[k]
                #phi <- phi / sqrt(sum(phi_new^2))
                h_old <- h_new
                if(h_new < h_min){
                    h_min <- h_new
                    phi_min <- phi / sqrt(sum(phi^2))
                }
            }
        }
        T <- T*0.9
    }
})
system.time({
    for(i in 1:(3000)){
        rand <- rnorm(36, mean = 0, sd = 0.2)*2
        rand2 <- runif(36,1)
        for(k in 1:36){
            h_new <- (phi %*% H %*% phi)/sum(phi^2)
        }
        T <- T*0.9
    }
})
print(h(H,phi_min))
print(phi_min)
eigen(H)$values
