# Toy example to debug

library(mvtnorm)

vecchia <- F
pmx <- T
one_layer <- F

if(one_layer) {source("../dgp.hm/R/logl_cov_1L.R")} else {source("../dgp.hm/R/logl_cov.R")}
source("../dgp.hm/R/bohman.R")
source("../dgp.hm/R/matrix.Moore.Penrose.R")
source("../dgp.hm/R/plot_fns.R") #plot.krig, plot.true, plot.warp
source("../dgp.hm/R/predict.R")
source("../dgp.hm/R/trim.R")
source("../dgp.hm/R/vecchia.R")

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

# Plots of observations, average, and true mean
# matplot(x, t(Y), col="gray", type="l")
# lines(x, ybar, col="red")
# lines(x, S_true, col="blue")

fit <- fit_two_layer_SW(x, ybar, nmcmc = 1000, Sigma_hat = Sigma_hat, pmx = pmx, vecchia = vecchia)
plot(fit)

fit <- trim_SW(fit, burn = 200, thin = 4)
plot(fit)

v <- fit$v
ytrue <- S_true
plot.true(fit, nrun = r)
plot.true.tau(fit, nrun = r, unif_tau = .1)
plot.warp(fit)
lines(x,w)
# plot.true.tau(fit, nrun = r, unif_tau = 1)
