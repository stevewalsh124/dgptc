# hetGP attempt to estimate sample covariance

library(mvtnorm)
library(MASS)
library(hetGP)

# this is from 1D_toy
r <- 16 #number of runs
n <- 100 #number of points per run
x <- seq(0,1,length.out=n)
w <- x^2

# true mean and covariance
S_true <- sin(4*pi*w)
Sigma_true <- diag(w^3 + .Machine$double.eps)

# generate r runs, each of length n
Y <- matrix(NA, r, n)
for(i in 1:r) Y[i,] <- S_true + rmvnorm(n = 1, sigma = Sigma_true)

# estimated mean and covariance
ybar <- colMeans(Y)
Sigma_hat <- cov(Y)/r

X <- matrix(rep(x, r, ncol=1))
Z <- c(t(Y))
# nvar <- 1
plot(X, Z, ylab = 'acceleration', xlab = "time", col="gray")


## Model fitting
model <- mleHetGP(X = X, Z = Z, #lower = rep(0.1, nvar), upper = rep(50, nvar),
                  covtype = "Gaussian")

## Display averaged observations
points(model$X0, model$Z0, pch = 20)

## A quick view of the fit                  
summary(model)

## Create a prediction grid and obtain predictions
xgrid <- matrix(x, ncol = 1)
predictions <- predict(x = xgrid, object =  model)

## Display mean predictive surface
lines(xgrid, predictions$mean, col = 'red', lwd = 2)
## Display 95% confidence intervals
lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
## Display 95% prediction intervals
lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
      col = 3, lty = 2)
lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
      col = 3, lty = 2)
## Display true function
lines(x, S_true, col="blue")

## Heteroscedastic nugget estimation. vs truth
plot(x, diag(Sigma_true), col="blue", type="l")
lines(x, model$Lambda)
