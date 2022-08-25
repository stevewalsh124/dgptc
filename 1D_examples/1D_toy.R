# Toy example to debug

library(mvtnorm)

vecchia <- T
pmx <- F
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
S_true <- sin(2*pi*w)
Sigma_true <- diag(w^3)

# generate r runs and find the sample average
Y <- matrix(NA, r, n)
for(i in 1:r) Y[i,] <- S_true + rmvnorm(n = 1, sigma = Sigma_true)
ybar <- colMeans(Y)

# Plots of observations, average, and true mean
# matplot(x, t(Y), col="gray", type="l")
# lines(x, ybar, col="red")
# lines(x, S_true, col="blue")

fit <- fit_two_layer_SW(x, ybar, nmcmc = 10000, Sigma_hat = cov(Y), pmx = pmx, vecchia = vecchia)
plot(fit)

fit <- trim_SW(fit, burn = 250, thin = 1)
plot(fit)
