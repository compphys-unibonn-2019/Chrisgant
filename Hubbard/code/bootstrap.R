N <- 50000 ## anzahl datenpunkte

x <- rnorm(n=N, mean=5, sd=0.1) ## x sollten deine daten sein, ist hier normalverzeiltes set

y <- sample.int(N, replace=TRUE) 

xstar <- x[y]                   ##xstar ist ein bootstrap sample, enthält die alle daten zufallspermutiert

R <- 200                        ##wenn das nicht zu lange dauert, sollte man das noch erhöhen
y <- array(sample.int(n=N, size=R*N, replace=TRUE), dim=c(R,N)) ##array aus R bootstrap samples mit je N einträgen
xstar  <- array(x[y], dim=c(R, N))

xbarstar <- apply(X=xstar, MARGIN=1, FUN=mean)          ##bestimme mittelwert für einzelne bootstrap samples
meanx <- mean(x)                                        ##mittelwert
deltax <- sd(xbarstar)                                  ##standard abweichung von mittelwert

bootstrap <- function(data = rnorm(n=1000, mean=5, sd=0.1), R=2000){
    N <- length(data)
    y <- sample.int(N, replace=TRUE) 
    #xstar <- data[y]
    y <- array(sample.int(n=N, size=R*N, replace=TRUE), dim=c(R,N))
    xstar  <- array(data[y], dim=c(R, N))
    xbarstar <- apply(X=xstar, MARGIN=1, FUN=mean)
    meanx <- mean(x)
    deltax <- sd(xbarstar)
    output <- list(mean=meanx, sd=deltax)
}

#qqnorm(xbarstar)
hist(xbarstar)
print(system.time({print(bootstrap(x,5000)$sd)}))