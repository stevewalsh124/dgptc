# Toy example to debug

library(mvtnorm)

vecchia <- F
pmx <- T
one_layer <- F

taper_cov <- T
tau_b <- 0.1
nmcmc <- 101000
nburn <- 1000
kth <- 5

if(one_layer) {source("../dgp.hm/R/logl_cov_1L.R")} else {source("../dgp.hm/R/logl_cov.R")}
source("../dgp.hm/R/bohman.R")
source("../dgp.hm/R/matrix.Moore.Penrose.R")
source("../dgp.hm/R/plot_fns.R") #plot.krig, plot.true, plot.warp
source("../dgp.hm/R/predict.R")
source("../dgp.hm/R/trim.R")
source("../dgp.hm/R/vecchia.R")

r <- 100 #number of runs
n <- 100 #number of points per run
x <- seq(0,1,length.out=n)
w <- x^3

pdf(paste0("pdf/1D_toy_",r,"_",n,"_",nmcmc,".pdf"))

# true mean and covariance
S_true <- sin(4*pi*w)
Sigma_true <- diag(w^3 + .Machine$double.eps)

# generate r runs, each of length n
Y <- matrix(NA, r, n)
for(i in 1:r) Y[i,] <- S_true + rmvnorm(n = 1, sigma = Sigma_true)

# estimated mean and covariance
ybar <- colMeans(Y)
Sigma_hat <- cov(Y)/r

# Plots of observations, average, and true mean
# matplot(x, t(Y), col="gray", type="l")
# lines(x, ybar, col="red")
# lines(x, S_true, col="blue")

# Fit and plots with toggle for tapering the covariance matrix
if(one_layer){
  if(taper_cov){
    fit <- fit_two_layer_SW(x, ybar, nmcmc = nmcmc, 
                            Sigma_hat = Sigma_hat * bohman(sqrt(plgp::distance(x)), tau_b))
  } else {
    fit <- fit_two_layer_SW(x, ybar, nmcmc = nmcmc, Sigma_hat = Sigma_hat)
  }
} else {
  if(taper_cov){
    fit <- fit_two_layer_SW(x, ybar, nmcmc = nmcmc, 
                            Sigma_hat = Sigma_hat * bohman(sqrt(plgp::distance(x)), tau_b), pmx = pmx, vecchia = vecchia)
  } else {
    fit <- fit_two_layer_SW(x, ybar, nmcmc = nmcmc, Sigma_hat = Sigma_hat, pmx = pmx, vecchia = vecchia)
  }
  
}

# trace plots before and after removing burn-in
plot(fit)
fit <- trim_SW(fit, burn = nburn, thin = kth)
plot(fit)

# look at UQ, via sims of the truth given the data
v <- fit$v
ytrue <- S_true
plot.true(fit, nrun = r)

# plot the hidden layer
if(!one_layer){
  plot.warp(fit)
  lines(x,w)
}

# study the (in)consistency/microergodicity of individual lengthscales vs their product
par(mfrow=c(1,3))
plot(fit$theta_w, type="l")
plot(fit$theta_y, type="l")
plot(fit$theta_y * fit$theta_w, type="l")

dev.off()
