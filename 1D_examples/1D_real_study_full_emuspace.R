# real data study

tic <- proc.time()[3]

saveImage <- T
PDF <- T

vecchia <- F
pmx <- F
one_layer <- F
force_id_warp <- F
load("w0_from_mte1_50k.rda")

# Smooth the precision info so it's not a step function for each data type (pt, lo, hi)?
smooth_precs <- T
if(smooth_precs) k_sm <- 10 #rolling mean uses k numbers

if(one_layer) {source("../dgp.hm/R/logl_cov_1L.R")} else {source("../dgp.hm/R/logl_cov.R")}
source("../dgp.hm/R/bohman.R")
source("../dgp.hm/R/matrix.Moore.Penrose.R")
source("../dgp.hm/R/plot_fns.R") #plot.krig, plot.true, plot.warp
source("../dgp.hm/R/predict.R")
source("../dgp.hm/R/trim.R")
source("../dgp.hm/R/vecchia.R")

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)
library(zoo) #rollmean
library(Matrix)
# library(latex2exp) #TeX in axis label

# Taper the covariance matrix before the MCMC fit?
taper_cov <- F

# Do a kriging step?
krig <- F

# Use hi res in Ybar calculation?
use_hi <- T

# Model the correlated errors with a covariance function?
cf_errors <- T
if(cf_errors){
  err_cov <- "matern"#"exp2"#
  loess_span <- 0.15
  err_v   <- paste0("est", loess_span)#ifelse(err_cov == "matern", 2.5, 999)
  err_g   <- NULL#sqrt(.Machine$double.eps)#
  err_g_msg <- ifelse(is.null(err_g),"estg","fixg")
}

ncores <- 2
tolpower <- -10

cov_fn <- "matern"#"exp2"#

if(taper_cov) tau_b <- .2
nrun <- 16
nmcmc <- 15000
nburn <- 10000
kth <- 4

bte <- 3 # cols 3-18 are low res

# Model 1, choose from 1-111 for training set, or c(0, 112:116) for test set
mte <- 1

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

if(PDF) pdf(paste0("pdf/realstudydgpact_",nmcmc,"_",nrun,"_",
                   if(cf_errors){paste0("_cfe",err_v,err_g_msg)},
                   if(taper_cov){paste0("tpr",tau_b)},
                   if(use_hi){"_uh"},if(pmx){"_pmx"},if(force_id_warp){"_fiw"},if(vecchia){"_vec"},
                   if(one_layer){"_1L"},"_",cov_fn,tolpower,"model",mte,"_",k_sm,".pdf"))

step <- 499

# first col is k, 2nd is linear pert theory, 3:18 is lowres, 19 is hires
if(mte %in% 0:111){
  pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                           if(mte<100){"0"},if(mte<10){"0"},mte,"_test.dat"))
} else {
  if(mte==112) pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",
                                        step,"/pk_E001_test.dat"))
  if(mte==113) pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",
                                        step,"/pk_E002_test.dat"))
  if(mte==114) pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",
                                        step,"/pk_E003_test.dat"))
  if(mte==115) pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",
                                        step,"/pk_E009_test.dat"))
  if(mte==116) pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",
                                        step,"/pk_E010_test.dat"))
}


n <- nrow(pk2)
k_sub <- 1:nrow(pk2) #index_list$lowres.ix

precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0) * nrun
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)

if(smooth_precs){
  Lam_lo <- diag(rollmean(precs_lo, k = k_sm, fill = "extend"))
  Lam_hi <- diag(rollmean(precs_hi, k = k_sm, fill = "extend"))
  Lam_pt <- diag(rollmean(precs_pt, k = k_sm, fill = "extend"))
} else {
  Lam_lo <- diag(precs_lo)
  Lam_hi <- diag(precs_hi)
  Lam_pt <- diag(precs_pt)
}

# subtract a linear model
# temp_avg <- rowMeans(log10(pk2[k_sub,3:18]))
# temp_lm <- lm(temp_avg ~ log10(pk2[k_sub,1]))$fitted.values
load("rda/avg_of_wt_avg_10.rda")
temp_lm <- rep(0,length(k_sub)) #log10(orig_loess*k^1.5/(2*pi^2))[k_sub]
#avg_loess[k_sub]
matplot(log10(k),log10(pk2[,3:18]))
matplot(log10(k)[k_sub],log10(pk2[k_sub,3:18]))

Y <- t(apply(t(pk2[k_sub,3:18]),1,function(x){log10(x*(k^1.5)[k_sub]/(2*pi^2))}))
Y <- ifelse(is.infinite(Y), 0, Y)

# wavenumber is X, a particular lowres run in Y
y_lo_avg <- colMeans(Y - temp_lm)
kl <- k[k_sub]
x <- log10(kl)
y <- log10(pk2[k_sub, bte]*kl^1.5/(2*pi^2)) - temp_lm
y_pt <- log10(pk2[k_sub, 2]*kl^1.5/(2*pi^2)) - temp_lm
y_pt <- ifelse(is.infinite(y_pt), 0, y_pt)
y_lra <- y_lo_avg#log10(y_lo_avg*k^1.5/(2*pi^2))[k_sub] - temp_lm
y_hi <- log10(pk2[k_sub, 19]*(k^1.5)[k_sub]/(2*pi^2)) - temp_lm
y_hi <- ifelse(is.infinite(y_hi), 0, y_hi)

# get a weighted average across low, high and pert
Lam_z <- Lam_pt+Lam_lo+Lam_hi
precs <- diag(Lam_z)
Lam_zi <- solve(Lam_z)
mu_z <- Lam_zi %*% (Lam_pt %*% y_pt + Lam_lo %*% t(t(y_lo_avg)) + Lam_hi %*% y_hi)

matplot(log10(k)[k_sub], t(Y), col="gray",type="l",lty=1, ylim = range(y_pt, y_lra, y_hi, Y))
lines(x, y_pt, col="green")
lines(x, y_lra, col="blue")
lines(x, y_hi, col="red")
lines(x[index_list$pert.ix], rep(-1, length(index_list$pert.ix)), col="green", lty=3)
lines(x[index_list$lowres.ix], rep(-1.1, length(index_list$lowres.ix)), col="blue",lty=3)
lines(x[index_list$highres.ix], rep(-1.2, length(index_list$highres.ix)), col = "red",lty=3)
lines(x, mu_z, col="orange",lty=2, lwd=4)

# hi res precs on the low res index (for comparison)
precs_hi <- prec_highres[k_sub]
if(use_hi){
  # Weighted average of low res runs and 
  hi_wt <- unique(prec_highres/prec_lowres)[1]
  y_avg <- mu_z#(hi_wt*y_hi+nrun*y_lra)/(hi_wt+nrun)
  nrunn <- nrun + hi_wt
} else {
  y_avg <- mu_z#y_lra
  nrunn <- nrun
}

# wavenumber is X, a particular lowres run in Y
mean_sz <- mean(y_avg)
sd_sz <- sd(y_avg)
x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean_sz)/sd_sz
y_avg <- (y_avg - mean_sz)/sd_sz
y_lra <- (y_lra - mean_sz)/sd_sz
y_hi <- (y_hi - mean_sz)/sd_sz
for (i in 1:nrow(Y)) { Y[i,] <- (Y[i,] - mean_sz)/sd_sz  }

D <- plgp:::distance(x)

# sd_sz is based on real data (sd for standardizing)
precc <- precs*sd_sz^2
sdd <- sqrt(1/precc)
A <- diag(sdd)

if(krig){
  # get predictions for xx and have corresponding prec info
  xx <- setdiff(seq(0,1,by=.01), x)
  lmfit <- lm(log10(precs) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
  betahat <- coef(lmfit)
  precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))
}

if(use_hi){
  n <- ncol(Y)
  r <- nrow(Y)
  # Ybar <- (r*y_lra+hi_wt*y_hi)/(r+hi_wt); same as y_avg above
  
  sum_YYt <- matrix(0, n, n)
  for (i in 1:r) { sum_YYt <- sum_YYt + (Y[i,]-y_avg)%*%t(Y[i,]-y_avg) }
  sum_YYt <- sum_YYt + hi_wt*(y_hi-y_avg)%*%t(y_hi-y_avg)
  # image.plot(cov(Y), zlim = range(c(cov(Y),1/(r+hi_wt-1)*sum_YYt)))
  # image.plot(1/(r+hi_wt-1)*sum_YYt, zlim = range(c(cov(Y),1/(r+hi_wt-1)*sum_YYt)))  
  Sigma_hat <- 1/(r+hi_wt-1) * sum_YYt
} else {
  Sigma_hat <- cov(Y)
}

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
    loess_fit <- loess(y_avg ~ x, span = loess_span)
    avg_loess <- loess_fit$fitted
  } else {
    avg_loess <- y_avg
  }
  
  Y_zm <- matrix(NA, nrow(Y), ncol(Y))
  for (ii in 1:nrow(Y_zm)) Y_zm[ii,] <- Y[ii,] - avg_loess
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
  
  lo_ind <- index_list$lowres.ix
  hi_ind <- index_list$highres.ix
  pt_ind <- index_list$pert.ix
  hi_only <- hi_ind[which(!(hi_ind %in% lo_ind))]
  
  out <- optim(c(.1, -5, .1, 1), nl_matern, method="L-BFGS-B", lower=c(-5,-35,-5,0.2),
               upper=c(4,3,4,100), D=D[lo_ind, lo_ind], Y=Y_zm[,lo_ind], A=A[lo_ind,lo_ind])
  
  # out <- optim(c(log(theta_e_true), log(g_e_true), log(tau2_e_true), true_v), 
  #              nl_matern, method="L-BFGS-B", lower=c(-5,-15,-5,0.2),
  #              upper=c(5,5,5,40), D=D, Y=Y_sim, A=A)
  
  kappa_hat <- out$par[4]
  phi_hat <- exp(out$par[1])
  theta_hat <- phi_hat*sqrt(kappa_hat)
  g_hat <- exp(out$par[2])
  tau2_hat <- exp(out$par[3])
  
  print(rbind(c("theta_e", "g_e", "tau2_e", "smooth"),
              c(theta_hat, g_hat, tau2_hat, kappa_hat)))
  
  Matern_hat <- tau2_hat * (geoR::matern(sqrt(D[lo_ind,lo_ind]), phi = phi_hat, kappa = kappa_hat) + diag(g_hat,ncol(Y[,lo_ind])))
  Sigma_hatty <-  A[lo_ind,lo_ind] %*% Matern_hat %*% A[lo_ind,lo_ind]
  # doubt this is right... need to redo the A matrices and precision info vectors
  Sigma_hat <- as.matrix(bdiag(solve(A[pt_ind,pt_ind])^0.5 %*% diag(rep(1/10000, length(pt_ind))) %*% solve(A[pt_ind,pt_ind])^0.5,
                     Sigma_hatty,
                     solve(A[hi_only,hi_only])^0.5 %*% diag((precs_hi[hi_only]*sd_sz^2)^-1) %*% solve(A[hi_only,hi_only])^0.5)) 
}

if(taper_cov) Sigma_hat <- Sigma_hat * bohman(sqrt(D), tau=tau_b)

par(mfrow=c(1,2))
matplot(x,t(Y),type="l",col="gray", ylim = range(y_hi, Y, y_avg, y_pt))
lines(x,y_hi,col="blue",lwd=2)
lines(x,y_avg,col="red",lwd=2)
legend("bottomright", legend = c("hi","avg"), lty=1, col=c("blue","red"))
image.plot(Sigma_hat/nrunn, main = "input as sigma_hat")

####################
# run for sim data #
####################

if(force_id_warp){
  fitcov <- fit_two_layer_SW(x = x, y = c(y_avg), nmcmc = nmcmc, Sigma_hat = Sigma_hat/nrunn, cov = cov_fn, pmx = pmx, 
                             vecchia = vecchia, settings = list(alpha =list(theta_w=1000), beta=list(theta_w=.0001/1000)))
} else {
  fitcov <- fit_two_layer_SW(x = x, y = c(y_avg), nmcmc = nmcmc, w_0 = w0_from_mte1_50k, 
                             Sigma_hat = Sigma_hat/nrunn, cov = cov_fn, pmx = pmx, 
                             vecchia = vecchia)
}

fitcov <- trim_SW(fitcov, nburn, kth)
plot(fitcov)
if(krig) {
  fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=ncores, precs_pred = precs_pred)
  
  ## MASSIVE HACK
  fitcov$s2 <- ifelse(fitcov$s2 < 0, 0, fitcov$s2)
  fitcov$s2_smooth <- ifelse(fitcov$s2_smooth < 0, 0, fitcov$s2_smooth)
  
  par(mfrow=c(1,1))
  plot.krig(fitcov)
}

v <- fitcov$v

par(mfrow=c(1,1))
fitcov <- est.true(fitcov)
plot.true(fitcov)
if(mte %in% 1:111){
  cosmicEmu <- read.csv(paste0("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/orig_111/EMU",
                               mte-1,".txt"),sep="", header = F)
} else {
  if(mte==0) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU5.txt",
                                   sep="", header = F)
  if(mte==112) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU0.txt",
                                     sep="", header = F)
  if(mte==113) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU1.txt",
                                     sep="", header = F)
  if(mte==114) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU2.txt",
                                     sep="", header = F)
  if(mte==115) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU3.txt",
                                     sep="", header = F)
  if(mte==116) cosmicEmu <- read.csv("~/CosmicEmu/2022-Mira-Titan-IV/P_tot/test_6/EMU4.txt",
                                     sep="", header = F)

}

cosmscrP <- log10(cosmicEmu[,2]*k^1.5/(2*pi^2))
cosmscrPsz <- (cosmscrP - mean_sz)/sd_sz
lines(fitcov$x, cosmscrPsz, col="green", lwd=1)

# Model 106 is hell
cb_cols   <- c("#999999", "#009E73", "#56B4E9", "#D55E00",
               "#E69F00", "#F0E442", "#0072B2", "#CC79A7")

# plot log10 space
# png(paste0("png/model",mte,"_emuspace.png"), width=4000, height = 2400, res=400)
par(mar=c(4,4.5,1,1), mfrow=c(1,1))
plot(fitcov$x, fitcov$y, type="n", #xlim = log10(c(.04,.35)),
     ylim = range(Y),
     xlab=expression(log[10](k)),
     ylab='script P')#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
for (i in 1:16) lines(fitcov$x, Y[i,], col="gray", lwd=1)
lines(fitcov$x, colMeans(Y), col=cb_cols[2], lwd=2)
lines(fitcov$x, y_hi, col=cb_cols[3], lwd=2)
# lines(fitcov$x, fitcov$Ms[,1], col="red", lwd=2)
lines(fitcov$x, fitcov$ubb, col=cb_cols[4], lwd=2, lty=3)
lines(fitcov$x, fitcov$lbb, col=cb_cols[4], lwd=2, lty=3)
lines(fitcov$x, cosmscrPsz, col=cb_cols[8], lwd=2, lty=2)
legend("bottomright",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
       col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
# dev.off()

m <- fitcov$m

# png(paste0("png/model",mte,"_emuspace_rmAvg.png"), width=4000, height = 2400, res=400)
par(mar=c(4,4.5,1,1), mfrow=c(1,1))
plot(log10(kl), fitcov$y - m, type="n", #xlim = log10(c(.04,.35)),
     ylim = c(-.33,.33),
     xlab=expression(log[10](k)),
     ylab='script P')#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
for (i in 1:16) lines(log10(kl), Y[i,] - m, col="gray", lwd=1)
lines(log10(kl), colMeans(Y) - m, col=cb_cols[2], lwd=2)
lines(log10(kl), y_hi - m, col=cb_cols[3], lwd=2)
# lines(log10(kl), fitcov$Ms[,1], col="red", lwd=2)
lines(log10(kl), fitcov$ubb - m, col=cb_cols[4], lwd=2, lty=3)
lines(log10(kl), fitcov$lbb - m, col=cb_cols[4], lwd=2, lty=3)
lines(log10(kl), cosmscrPsz - m, col=cb_cols[8], lwd=2, lty=2)
legend("topleft",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
       col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
# dev.off()


if(!one_layer) plot.warp(fitcov)

if(PDF) dev.off()

toc <- proc.time()[3]

(timing <- toc - tic)

if(saveImage) save.image(file = paste0("/projects/precipit/1D_real_study/rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                                       if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                                       if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))
