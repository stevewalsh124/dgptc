############################################
# Simulating deep GPs, naively
# Steve Walsh
# January 28 2022
############################################

# Code based on that of M.A.R. Ferreira, 2016.

library(fields)

################################
# Gaussian covariance function #
################################

covfunc.Gaussian <- function(t,phi,sigma2) {sigma2 * exp(-t^2/phi^2)}

matrix.sqrt <- function(H)
{
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-10)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  

  H.sqrt
}

# With nugget effect
plot.gaus <- function(t, phi=1, sigma2=1, tau2=0){
  covmatrix <- diag(tau2,nrow(t)) + covfunc.Gaussian(t,phi,sigma2)
  x <- matrix.sqrt(covmatrix) %*% rnorm(nrow(s))
  z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
  image.plot(z, main = bquote(phi == .(phi) ~ ",  " ~ sigma^2 == 
                                .(sigma2) ~ ",  " ~ tau^2 == .(tau2)),
             cex.main=1.5)
}


####################
# Plot Simulations #
####################

par(mfrow=c(2,3))

# Currently the same for each layer; can change for individual layers
phi <- 1/3
sigma2 <- 1
tau2 <- .Machine$double.eps^0.5

xaxis <- seq(0,1,0.04)
yaxis <- seq(0,1,0.04)

# original lat lon coords
s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
# plot(s)

# first layer
t <- as.matrix(dist(s))
covmatrix <- diag(tau2,nrow(s)) + covfunc.Gaussian(t,phi,sigma2)
x <- t(chol(covmatrix)) %*% rnorm(nrow(s))

# second layer
t2 <- as.matrix(dist(x))
covmatrix2 <- diag(tau2,nrow(s)) + covfunc.Gaussian(t2,phi,sigma2)
w <- t(chol(covmatrix2)) %*% rnorm(nrow(x))

# # third layer
# t3 <- as.matrix(dist(w))
# covmatrix3 <- diag(tau2,nrow(s)) + covfunc.Gaussian(t3,phi,sigma2)
# z <- t(chol(covmatrix2)) %*% rnorm(nrow(w))

# output; change "matrix(w" to "matrix(z" if using third layer
y <- matrix(w,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
nug <- rnorm(prod(dim(y)), sd = sqrt(tau2))
image.plot(y + matrix(nug, dim(y)[1], dim(y)[2]), 
           main = bquote(phi == .(phi) ~ ",  " ~ sigma^2 == .(sigma2) ~ ",  "
                         ~ tau^2 == .(round(tau2,4))), cex.main=1.5)