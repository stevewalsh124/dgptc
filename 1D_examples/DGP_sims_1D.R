############################################
# Simulating deep GPs, naively
# Steve Walsh
# January 28 2022
############################################

# Based on Bobby Gramacy's "Surrogates", Chapter 5

n <- 500
X <- matrix(seq(0, 10, length=n), ncol=1)

library(plgp) #plgp::distance
library(mvtnorm)

#squared distance
D <- distance(X) 

# lengthscale = 1
# marginal variance = 1
eps <- sqrt(.Machine$double.eps) 
Sigma <- exp(-D) + diag(eps, n)

# first layer
W <- c(rmvnorm(1, mean=X, sigma=Sigma))

# second layer
D2 <- distance(c(W))
Sigma_2 <- exp(-D2) + diag(eps, n)
Y <- rmvnorm(1, sigma=Sigma_2) #mean=W,

#Plot each combination of layers/output
par(mfrow=c(2,2))
plot(X, W, type="l"); abline(0, 1, col=2, lty=2)
plot(W, Y, type = "l"); abline(h=0, col=2, lty=2)
plot(X, Y, type="l"); abline(h=0, col=2, lty=2)
