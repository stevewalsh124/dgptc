# sim study

tic <- proc.time()[3]

saveImage <- T
PDF <- T

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
# source("get_matern.R")

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
model_diag <- T
var_adj <- 40000 #23.55 is MLE from c_hat.R

# Model the correlated errors with a covariance function?
cf_errors <- F
if(cf_errors){
  err_cov <- "matern"#"exp2"#
  loess_span <- 0.15
  err_v   <- paste0("est", loess_span)#ifelse(err_cov == "matern", 2.5, 999)
  err_g   <- NULL#sqrt(.Machine$double.eps)#
  err_g_msg <- ifelse(is.null(err_g),"estg","fixg")
}

# Use the true covariance matrix for sigma hat?
use_true_cov <- F

# Use both the true cov for Sigma hat as well as Sigma_S
use_both_true_covs <- F
if(use_both_true_covs) use_true_cov <- T

# Taper the covariance matrix before the MCMC loop?
taper_cov <- F
tau_b <- .2

# Don't do this; use est.true* instead
krig <- F

seed <- 481

cov_fn <- "matern"#"exp2"#

# Number of computer experiment runs to simulate
nrun <- 16

# MCMC info, including thinning every kth
nmcmc <- 20000
nburn <- 10000
kth <- 4

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

if(PDF) pdf(paste0("pdf/simstudydgpact_",nmcmc,"_",nrun,
           if(true_diag){"_TD"}, if(model_diag){paste0("_MD", var_adj)},
           if(cf_errors){paste0("_cfe",err_v,err_g_msg)},
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

# lengthscale = theta_w_true
# marginal variance = 1
if(seed <= 100){
  theta_w_true <- 1
  Sigma_W_true <- exp(-D/theta_w_true) + diag(eps, length(x))
  # first layer
  # load("rda/warp_true.rda")
  w <- c(rmvnorm(1, mean=x, sigma=Sigma_W_true))#warp_true#
  
  # second layer
  Dw <- plgp:::distance(c(w))
  theta_y_true <- 0.05
  tau2_true <- 1
  g_true <- eps
  Sigma_S_true <- tau2_true * exp(-Dw/theta_y_true) + diag(tau2_true * g_true, length(w))
  ytrue_dgp <- rmvnorm(1, sigma=Sigma_S_true) #mean=W,
  ytrue_dgp <- (ytrue_dgp - mean(ytrue_dgp))/sd(ytrue_dgp)
  
  #Plot each combination of layers/output
  par(mfrow=c(1,3))
  plot(x, w, type="l"); abline(0, 1, col=2, lty=2)
  plot(w, ytrue_dgp, type = "l"); abline(h=0, col=2, lty=2)
  plot(x, ytrue_dgp, type="l"); abline(h=0, col=2, lty=2)
}
# lmf <- lm(y_avg ~ x)$fitted.values
# ytrue <- y_avg - lmf
# ytrue <- (ytrue - mean(ytrue))/sd(ytrue)
# plot(ytrue, type="l", main = "true y")

# generate 16 low res runs with iid noise (to start)
# Gaussian cov fn with theta=1, tau2=1

# sd_sz is based on real data (sd for standardizing)
if(seed <= 400) sd_sz <- 0.009034001 else sd_sz <- 1
precc <- precs*sd_sz^2
sdd <- sqrt(1/precc)

if(true_diag){
  Sigma_e_true <- diag(1/precc)
} else {
  if(seed <= 100) {
    theta_e_true <- 0.01 
    tau2_e_true <- 2
    load("rda/yavg_for_sims.rda")
    ytrue <- c(ytrue_dgp) #+ y_avg_true
    }else{
      theta_e_true <- 0.0005
      tau2_e_true <- .25
      if(seed > 100 & seed <= 200) {
        load("rda/loess_scriptPs_0.1.rda")
        ytrue <- (loess_scriptPs[seed-100,] - mean(loess_scriptPs[seed-100,]))/sd(loess_scriptPs[seed-100,])
      } else { # seed > 200
        if(seed > 200 & seed <= 311){
          # CosmicEmu makes .txt files 0-110 rather than 1-100, so seeds 201-311 go to EMU 0-110
          ytrue_orig <- read.csv(paste0("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/orig_111/EMU",
                                        seed-201,".txt"), sep = "", header = F)
          ytrue_un <- log10(ytrue_orig[,1]^1.5*ytrue_orig[,2]/(2*pi^2))[index_list$lowres.ix]
          ytrue_sd <- sd(ytrue_un); ytrue_mean <- mean(ytrue_un)
          ytrue <- (ytrue_un - ytrue_mean) / ytrue_sd
        } else { #seed > 311
          if(seed > 311 & seed < 401) stop("forbidden seed")
          if(seed > 511) stop("forbidden seed")
          # CosmicEmu makes .txt files 0-110 rather than 1-100, so seeds 201-311 go to EMU 0-110
          ytrue_orig <- read.csv(paste0("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/orig_111/EMU",
                                        seed-401,".txt"), sep = "", header = F)
          ytrue_un <- log10(ytrue_orig[,1]^1.5*ytrue_orig[,2]/(2*pi^2))[index_list$lowres.ix]
          ytrue_sd <- sd(ytrue_un); ytrue_mean <- mean(ytrue_un)
          ytrue <- (ytrue_un - ytrue_mean) / ytrue_sd
        }
      }
    }
  if(seed <= 311){
    g_e_true <- eps
    Sig_e_temp <- tau2_e_true * exp(-D/theta_e_true) + diag(tau2_e_true * g_e_true, length(x)) 
  } else {
    g_e_true <- 2e-12
    tau2_e_true <- 2 / ytrue_sd^2
    smooth_true <- 4
    theta_e_true <- 0.02/sqrt(smooth_true)
    Sig_e_temp <- tau2_e_true * (geoR::matern(sqrt(D),phi = theta_e_true,kappa = smooth_true) + diag(g_e_true,length(x)))
  }
  A <- diag(sdd)
  Sigma_e_true <- A %*% Sig_e_temp %*% A
}

Y_sim <- mvtnorm::rmvnorm(n = nrun, mean = ytrue, sigma = Sigma_e_true, method = "eigen")

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
  Sigma_e_true_pred <- diag(sqrt(1/precs_pred)/sd_sz)
} else {
  Sig_e_temp_p <- exp(-plgp:::distance(xx/.05)) + diag(sqrt(.Machine$double.eps),length(xx)) 
  Ap <- diag(sqrt(1/precs_pred)/sd_sz)
  Sigma_e_true_pred <- Ap %*% Sig_e_temp_p %*% Ap
}

# get sigma hat
if(use_true_cov) {
  Sigma_hat <- Sigma_e_true
} else {
  if(model_diag){
    Sigma_hat <- var_adj*diag(1/precc)
  } else {
    if(cf_errors){
      # # logl_cov* files use the same names (eg: logl_SW, fit_two_layer_SW)
      # if(!one_layer) source("../dgp.hm/R/logl_cov_1L.R")
      # varvec <- 1/precc
      # Sigma_hat_ho <- get_matern(x, Y_sim - colMeans(Y_sim), nmcmc = err_mcmc, nburn = err_burn, 
      #                         cov = err_cov, v = err_v, true_g = err_g, varvec = varvec)
      # Sigma_hat <- diag(sqrt(varvec)) %*% Sigma_hat_ho %*% diag(sqrt(varvec))
      # if(!one_layer) source("../dgp.hm/R/logl_cov.R")
      # D <- plgp::distance(x)
      if(loess_span != 0){
        loess_fit <- loess(colMeans(Y_sim) ~ x, span = loess_span)
        avg_loess <- loess_fit$fitted
      } else {
        avg_loess <- colMeans(Y_sim)
      }
      
      Y_sim_zm <- matrix(NA, nrow(Y_sim), ncol(Y_sim))
      for (ii in 1:nrow(Y_sim_zm)) Y_sim_zm[ii,] <- Y_sim[ii,] - avg_loess
      Y <- Y_sim_zm
      counter <- 0
      
      nl_matern <- function(par, D, Y, A) {
        theta <- exp(par[1])
        g <- exp(par[2])
        tau2 <- exp(par[3])
        kappa <- par[4]
        if(counter %% 50 == 0) print(c(theta,g,tau2,kappa))
        n <- ncol(Y)
        K <- tau2 * (geoR::matern(sqrt(D),phi = theta,kappa = kappa) + diag(g,n))
        
        Ki <- solve(K)
        ldetK <- determinant(K, logarithm=TRUE)$modulus
        
        ll <- 0
        for (ii in 1:nrow(Y)) {
          yit <- solve(A) %*% Y[ii,]
          ll <- ll - (1/2)*(t(yit) %*% Ki %*% yit) - (1/2)*ldetK 
        }
        counter <<- counter + 1
        Kii <<- Ki
        return(-ll)
      }
      
      out <- optim(c(.1, -5, .1, 1), nl_matern, method="L-BFGS-B", lower=c(-5,-15,-5,0.2),
                   upper=c(10,1,10,100), D=D, Y=Y, A=A)
      
      kappa_hat <- out$par[4]
      phi_hat <- exp(out$par[1])
      theta_hat <- phi_hat*sqrt(kappa_hat)
      g_hat <- exp(out$par[2])
      tau2_hat <- exp(out$par[3])
      
      print(rbind(c("theta_e", "g_e", "tau2_e", "smooth"),
                  c(theta_hat, g_hat, tau2_hat, kappa_hat),
                  c(theta_e_true, g_e_true, tau2_e_true, "true_v")))
      
      Matern_hat <- tau2_hat * (geoR::matern(sqrt(D), phi = phi_hat, kappa = kappa_hat) + diag(g_hat,ncol(Y)))
      Sigma_hat <-  A %*% Matern_hat %*% A
    } else { 
      Sigma_hat <- cov(Y_sim)
    }
  }
}

if(taper_cov) Sigma_hat <- Sigma_hat * bohman(sqrt(plgp:::distance(x)), tau_b)

par(mfrow=c(1,2))
image.plot(Sigma_e_true, main = "Sigma_e_true", zlim = range(c(Sigma_e_true,Sigma_hat)))
image.plot(Sigma_hat, main = "input as sigma_hat", zlim = range(c(Sigma_e_true,Sigma_hat)))

####################
# run for sim data #
####################

if(one_layer){
  fitcov <- fit_two_layer_SW(x = x, y = colMeans(Y_sim), nmcmc = nmcmc, 
                             Sigma_hat = Sigma_hat/nrun, cov = cov_fn, vecchia = vecchia)
} else {
  fitcov <- fit_two_layer_SW(x = x, y = colMeans(Y_sim), nmcmc = nmcmc, 
                             Sigma_hat = Sigma_hat/nrun, cov = cov_fn,
                             pmx = pmx, vecchia = vecchia)
}

# plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
if(one_layer){
  par(mfrow=c(1,3))
  plot(fitcov$g, type="l", xlab="Iteration", ylab="g", main="Trace Plot of g")
  if(seed <= 100) {abline(h=g_true, col="blue"); abline(h=g_true/nrun, col="blue", lty=2)}
  plot(fitcov$theta_y, type="l", xlab="Iteration", ylab="theta_y", main="Trace Plot of theta_y")
  if(seed <= 100) abline(h=theta_y_true, col="blue")
  plot(fitcov$tau2, type="l", xlab="Iteration", ylab="tau2", main="Trace Plot of tau2")
  if(seed <= 100) abline(h=tau2_true, col="blue")
} else {
  par(mfrow=c(2,2))
  plot(fitcov$g, type="l", xlab="Iteration", ylab="g", main="Trace Plot of g")
  if(seed <= 100) {abline(h=g_true, col="blue"); abline(h=g_true/nrun, col="blue", lty=2)}
  plot(fitcov$theta_y, type="l", xlab="Iteration", ylab="theta_y", main="Trace Plot of theta_y")
  if(seed <= 100) abline(h=theta_y_true, col="blue")
  plot(fitcov$theta_w, type="l", xlab="Iteration", ylab="theta_w", main="Trace Plot of theta_w")
  if(seed <= 100) abline(h=theta_w_true, col="blue")
  plot(fitcov$tau2, type="l", xlab="Iteration", ylab="tau2", main="Trace Plot of tau2")
  if(seed <= 100) abline(h=tau2_true, col="blue")
  if(krig) fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=2, precs_pred = 1/diag(Sigma_e_true_pred))
}

v <- fitcov$v

par(mfrow=c(1,1))
if(krig) plot.krig(fitcov, Y = Y_sim)
if(use_both_true_covs){
  fitcov <- est.true.truecovs(fitcov, Y = Y_sim, Sigma_s = Sigma_S_true)
} else {
  fitcov <- est.true(fitcov, Y = Y_sim)
}
plot.true(fitcov, Y=Y_sim) # writes csvs of the emp coverage (for both tapering and no tapering)
# if(!taper_cov) plot.true.tau(fitcov, Y=Y_sim, unif_tau = tau_b) #no longer writes csvs

if(!one_layer){
  if(seed <= 100){
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
  }
  plot.warp(fitcov)
}

if(PDF) dev.off()

toc <- proc.time()[3]

timing <- toc - tic

if(saveImage){
  write.csv(timing, 
            file = paste0("csv/timing/timing_",
                          nmcmc,"_",nrun,
                          if(true_diag){"_TD"}, 
                          if(model_diag){paste0("_MD", var_adj)},
                          if(cf_errors){paste0("_cfe",err_v,err_g_msg)},
                          if(use_true_cov){"_UTC"}, if(use_both_true_covs){"_UBTC"},
                          if(taper_cov){paste0("_tpr",tau_b)}, 
                          if(pmx){"_pmx"}, if(vecchia){"_vec"},
                          if(one_layer){"_1L"},"_",cov_fn,"_",seed,".csv"))
  
  save.image(file = paste0("/projects/precipit/1D_sim_study/rda/1D_sim_",nmcmc,"_",nrun,
                           if(true_diag){"_TD"}, if(model_diag){paste0("_MD", var_adj)},
                           if(cf_errors){paste0("_cfe",err_v,err_g_msg)},
                           if(use_true_cov){"_UTC"}, if(use_both_true_covs){"_UBTC"},
                           if(taper_cov){paste0("_tpr",tau_b)}, if(pmx){"_pmx"}, if(vecchia){"_vec"},
                           if(one_layer){"_1L"},"_",cov_fn,"_",seed,".rda"))
}
