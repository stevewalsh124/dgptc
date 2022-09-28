# sim study

one_layer <- F
pmx <- F
vecchia <- F

if(one_layer) {source("../dgp.hm/R/logl_cov_1L.R")} else {source("../dgp.hm/R/logl_cov.R")}
source("../dgp.hm/R/matrix.Moore.Penrose.R")
source("../dgp.hm/R/plot_fns.R") #plot.krig, plot.true, plot.warp
source("../dgp.hm/R/bohman.R")
source("../dgp.hm/R/predict.R")
source("../dgp.hm/R/trim.R")
source("../dgp.hm/R/vecchia.R")

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
precs <- prec_lowres[index_list$lowres.ix]

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)

# Should the true error's cov mtx be diagonal or dense?
true_diag <- F

# Should the observations' noise be modeled as diagonal?
model_diag <- F
var_adj <- 25

# Use the true covariance matrix for sigma hat?
use_true_cov <- T

# Use both the true cov for Sigma hat as well as Sigma_S
use_both_true_covs <- T
if(use_both_true_covs) use_true_cov <- T

# Taper the covariance matrix before the MCMC loop?
taper_cov <- F
tau_b <- 1

# Don't do this; use est.true* instead
krig <- F

seed <- 1

cov_fn <- "matern"#"exp2"#

# Number of computer experiment runs to simulate
nrun <- 16

# MCMC info, including thinning every kth
nmcmc <- 12000
nburn <- 2000
kth <- 5

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

pdf(paste0("pdf/simstudydgpact_",nmcmc,"_",nrun,
           if(true_diag){"_TD"}, if(model_diag){paste0("_MD", var_adj)},
           if(use_true_cov){"_UTC"}, if(use_both_true_covs){"_UBTC"},
           if(taper_cov){paste0("_tpr",tau_b)}, if(pmx){"_pmx"}, if(vecchia){"_vec"},
           if(one_layer){"_1L"},"_",cov_fn,"_",seed,".pdf"))

# The step number, corresponding to a particular redshift
# steps <- c(163, 189, 247, 300, 347, 401, 453, 499)
step <- 499

# Model 1, choose from 000-111; exploring the 
# 8-dimensional cosmological parameter space
i <- 1

pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# hi res precs on the low res index (for comparison)
precs_hi <- prec_highres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
# bte <- 3 # band to evaluate; cols 3-18 are low res
# y <- log10(pk2[index_list$lowres.ix, bte])
# Y <- log10(pk2[index_list$lowres.ix, 3:18])
x <- (x - min(x))/(max(x)-min(x))
# y <- (y - mean(y))/sd(y)
# y_avg <- rowMeans(Y)
# y_avg <- (y_avg - mean(y_avg))/sd(y_avg)

set.seed(seed)
# computes squared distances
D <- plgp:::distance(x)

# lengthscale = 1
# marginal variance = 1
Sigma <- exp(-D) + diag(eps, length(x))

# first layer
set.seed(seed)
load("rda/warp_true.rda")
w <- warp_true#c(rmvnorm(1, mean=x, sigma=Sigma))#

# second layer
Dw <- plgp:::distance(c(w))
Sigma_s <- exp(-Dw/0.05) + diag(eps, length(w))
ytrue_dgp <- rmvnorm(1, sigma=Sigma_s) #mean=W,
ytrue_dgp <- (ytrue_dgp - mean(ytrue_dgp))/sd(ytrue_dgp)

#Plot each combination of layers/output
par(mfrow=c(1,3))
plot(x, w, type="l"); abline(0, 1, col=2, lty=2)
plot(w, ytrue_dgp, type = "l"); abline(h=0, col=2, lty=2)
plot(x, ytrue_dgp, type="l"); abline(h=0, col=2, lty=2)

# lmf <- lm(y_avg ~ x)$fitted.values
# ytrue <- y_avg - lmf
# ytrue <- (ytrue - mean(ytrue))/sd(ytrue)
# plot(ytrue, type="l", main = "true y")

# generate 16 low res runs with iid noise (to start)
# Gaussian cov fn with theta=1, tau2=1

# sd_sz is based on real data (sd for standardizing)
sd_sz <- 0.009034001
precc <- precs*sd_sz^2
sdd <- sqrt(1/precc)

if(true_diag){
  Cov_true <- diag(1/precc)
} else {
  Sig_mat <- exp(-plgp:::distance(x/.05)) + diag(sqrt(.Machine$double.eps),length(x)) 
  A <- diag(sdd)
  Cov_true <- A %*% Sig_mat %*% A
}

load("rda/yavg_for_sims.rda")
ytrue <- c(ytrue_dgp) #+ y_avg_true
Y_sim <- mvtnorm::rmvnorm(n = nrun, mean = ytrue, sigma = Cov_true, method = "eigen")
set.seed(as.numeric(Sys.time()))

par(mfrow=c(1,1))
matplot(x, t(Y_sim), type = "l", col="gray", main = "true vs estimated curve")
lines(x, ytrue, col="red", lwd=2)
lines(x,colMeans(Y_sim), lwd=2, lty=2)
legend("bottomright", legend = c("true", "sample avg"), col=c("red","black"), lty=1:2, lwd=2)


# get predictions for xx and have corresponding prec info
xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(precs) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))
if(true_diag){
  Sig_mat_p <- exp(-plgp:::distance(xx/.05)) + diag(sqrt(.Machine$double.eps),length(xx)) 
  Ap <- diag(2500/precs_pred^0.8)
  Cov_true_pred <- Ap %*% Sig_mat_p %*% Ap
} else {
  Cov_true_pred <- diag(5000/precs_pred^0.8)
}

# get sigma hat
if(use_true_cov) {
  Sigma_hat <- Cov_true
} else {
  if(model_diag){
    Sigma_hat <- var_adj*diag(1/precc)
  } else {
    Sigma_hat <- cov(Y_sim)
  }
}

if(taper_cov) Sigma_hat <- Sigma_hat * bohman(sqrt(plgp:::distance(x)), tau_b)

par(mfrow=c(1,2))
image.plot(Cov_true, main = "Cov_true", zlim = range(c(Cov_true,Sigma_hat)))
image.plot(Sigma_hat, main = "input as sigma_hat", zlim = range(c(Cov_true,Sigma_hat)))

####################
# run for sim data #
####################

fitcov <- fit_two_layer_SW(x = x, y = colMeans(Y_sim), nmcmc = nmcmc, 
                           Sigma_hat = Sigma_hat/nrun, cov = cov_fn,
                           pmx = pmx, vecchia = vecchia)
# plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
if(one_layer){
  plot(fitcov$g, type="l")
  plot(fitcov$theta_y, type = "l")
} else {
  plot(fitcov)
  fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=2, precs_pred = 1/diag(Cov_true_pred))
}

v <- fitcov$v

par(mfrow=c(1,1))
if(krig) plot.krig(fitcov, Y = Y_sim)
if(use_both_true_covs){
  fitcov <- est.true.truecovs(fitcov, Y = Y_sim, Sigma_s = Sigma_s)
} else {
  fitcov <- est.true(fitcov, Y = Y_sim)
}
plot.true(fitcov, Y=Y_sim) # writes csvs of the emp coverage (for both tapering and no tapering)
# if(!taper_cov) plot.true.tau(fitcov, Y=Y_sim, unif_tau = tau_b) #no longer writes csvs

if(!one_layer){
  # This modifies the true warping in a Procrustes-type transformation
  # So that it is rescaled, and rotated and with endpoints at 0 and 1
  wl = 1; wh = length(fitcov$x); ref.scale = 1

  # Second: undo scale
  # find the length of the deformed reference vector
  y.ref <- w[wl] - w[wh]
  scale.y <- sqrt(t(y.ref)%*%y.ref)
  # numerator comes from length of vector with points zz, zo
  scale.back =  c(ref.scale/scale.y)#1/sqrt(sum(y.translate[,2]^2))
  y.rot.scale = w * scale.back
  
  # Third: Translate back
  translation <- y.rot.scale[wl] - x[wl]
  y.rot.scale.tran <- y.rot.scale - translation
  warp_true <- y.rot.scale.tran
  
  # flip the graph if necessary in either direction
  # remove the reflections
  if(warp_true[wl]>warp_true[wh]) {warp_true <- (warp_true-warp_true[wl])/(warp_true[wh]-warp_true[wl])}
  plot.warp(fitcov)
}

dev.off()
