N <- 1000 ## anzahl datenpunkte

x <- rnorm(n=N, mean=5, sd=0.1) ## x sollten deine daten sein, ist hier normalverzeiltes set

y <- sample.int(N, replace=TRUE) 

xstar <- x[y]                   ##xstar ist ein bootstrap sample, enthält die alle daten zufallspermutiert

R <- 200                        ##wenn das nicht zu lange dauert, sollte man das noch erhöhen
y <- array(sample.int(n=N, size=R*N, replace=TRUE), dim=c(R,N)) ##array aus R bootstrap samples mit je N einträgen
xstar  <- array(x[y], dim=c(R, N))

xbarstar <- apply(X=xstar, MARGIN=1, FUN=mean)          ##bestimme mittelwert für einzelne bootstrap samples
meanx <- mean(x)                                        ##mittelwert
deltax <- sd(xbarstar)                                  ##standard abweichung von mittelwert


##qqnorm(xbarstar)
##hist(xbarstar)