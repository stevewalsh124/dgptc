# Combine Perturbation Theory, Low Res, and Hi Res runs
# See "combine_pert_lo_hi.pdf" for modeling details
# real data study

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

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)

# Taper the covariance matrix before the MCMC fit?
taper_cov <- F

# Do a kriging step?
krig <- F

# Use hi res in Ybar calculation?
use_hi <- T

ncores <- 2
tolpower <- -10

cov_fn <- "matern"#"exp2"#

tau_b <- 1
nrun <- 16
nmcmc <- 1000
nburn <- 200
kth <- 2

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

pdf(paste0("pdf/compare_pt_wts",nmcmc,"_",one_layer,if(pmx){"_pmx"},if(vecchia){"_vec"},".pdf"))
# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")

n <- length(k)

precs_lo <- ifelse(1:n %in% index_list$lowres.ix, prec_lowres, 0)
precs_hi <- ifelse(1:n %in% index_list$highres.ix, prec_highres, 0)
precs_pt <- ifelse(1:n %in% index_list$pert.ix, 10000, 0)

Lam_lo <- diag(precs_lo)
Lam_hi <- diag(precs_hi)
Lam_pt <- diag(precs_pt)

# step has to do with where/what time to which the simulation relates
step <- 499

# Model 1, choose from 000-111
i <- 1

# first col is k, 2nd is linear pert theory, 3:18 is lowres, 19 is hires
pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
y_pt <- pk2[,2]
y_lo_avg <- rowMeans(pk2[,3:18])
y_hi <- pk2[,19]
x_un <- log10(k)
x <- (x_un - min(x_un))/(max(x_un)-min(x_un))

plot(x, log10(y_lo_avg), type="l", 
     ylim = range(c(log10(y_pt[index_list$pert.ix]), log10(y_lo_avg[index_list$lowres.ix]), log10(y_hi[index_list$highres.ix]))))
# for (i in 3:18) lines(x, log10(pk2[,i]))
lines(x, log10(y_hi), col="red")
lines(x, log10(y_pt), col="orange")

# get a weighted average across low, high and pert
Lam_z <- Lam_pt+Lam_lo+Lam_hi
Lam_zi <- solve(Lam_z)
mu_z <- Lam_zi %*% (Lam_pt %*% y_pt + Lam_lo %*% t(t(y_lo_avg)) + Lam_hi %*% y_hi)
lines(x, log10(mu_z), col="blue",lty=2)

# pt_wt <- 0.00001
# precs_pt <- ifelse(1:n %in% index_list$pert.ix, pt_wt, 0)
# Lam_pt <- diag(precs_pt)
# Lam_z <- Lam_pt+Lam_lo+Lam_hi
# Lam_zi <- solve(Lam_z)
# mu_z <- Lam_zi %*% (Lam_pt %*% y_pt + Lam_lo %*% t(t(y_lo_avg)) + Lam_hi %*% y_hi)
# lines(x, log10(mu_z), col="green", lty=3)

# subtract a linear model
x2 <- x^2
temp_lm <- lm(log10(mu_z) ~ x + x2)$fitted.values
y_un <- c(log10(mu_z) - temp_lm)
y <- (y_un - mean(y_un))/sd(y_un)
lines(x, temp_lm, col="green")
legend("bottomleft",legend=c("pert","hi","lo","wt avg","lm"), col = c("orange","red","black","blue","green"), lty=c(1,1,1,2,1))
lines(x[index_list$pert.ix], rep(2.35, length(index_list$pert.ix)))
lines(x[index_list$lowres.ix], rep(2.5, length(index_list$lowres.ix)), lty=2)
lines(x[index_list$highres.ix], rep(2.65, length(index_list$highres.ix)), lty=3)

# recompute upon adjusting for scaling (by var(y_un)) of log10(P(k))
Lam_lo <- diag(precs_lo)*var(y_un)
Lam_hi <- diag(precs_hi)*var(y_un)
Lam_pt <- diag(precs_pt)*var(y_un)
Lam_z <- Lam_pt+Lam_lo+Lam_hi
Lam_zi <- solve(Lam_z)

# # wavenumber is X, a particular lowres run in Y
# x <- log10(k[index_list$lowres.ix])
# y <- log10(pk2[index_list$lowres.ix, bte]) - temp_lm
# y_lra <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]) - temp_lm)
# y_hi <- log10(pk2[index_list$lowres.ix, 19]) - temp_lm
# 
# # hi res precs on the low res index (for comparison)
# precs_hi <- prec_highres[index_list$lowres.ix]
# if(use_hi){
#   # Weighted average of low res runs and 
#   hi_wt <- unique(prec_highres/prec_lowres)[1]
#   y_avg <- (hi_wt*y_hi+nrun*y_lra)/(hi_wt+nrun)
#   nrunn <- nrun + hi_wt
# } else {
#   y_avg <- y_lra
#   nrunn <- nrun
# }
# 
# # wavenumber is X, a particular lowres run in Y
# mean_sz <- mean(y_avg)
# sd_sz <- sd(y_avg)
# x <- log10(k[index_list$lowres.ix])
# x <- (x - min(x))/(max(x)-min(x))
# y <- (y - mean_sz)/sd_sz
# y_avg <- (y_avg - mean_sz)/sd_sz
# y_lra <- (y_lra - mean_sz)/sd_sz
# y_hi <- (y_hi - mean_sz)/sd_sz
# for (i in 1:nrow(Y)) { Y[i,] <- (Y[i,] - mean_sz)/sd_sz  }
# 
# if(krig){
#   # get predictions for xx and have corresponding prec info
#   xx <- setdiff(seq(0,1,by=.01), x)
#   lmfit <- lm(log10(precs) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
#   betahat <- coef(lmfit)
#   precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))
# }
# 
# if(use_hi){
#   n <- ncol(Y)
#   r <- nrow(Y)
#   # Ybar <- (r*y_lra+hi_wt*y_hi)/(r+hi_wt); same as y_avg above
#   
#   sum_YYt <- matrix(0, n, n)
#   for (i in 1:r) { sum_YYt <- sum_YYt + (Y[i,]-y_avg)%*%t(Y[i,]-y_avg) }
#   sum_YYt <- sum_YYt + hi_wt*(y_hi-y_avg)%*%t(y_hi-y_avg)
#   # image.plot(cov(Y), zlim = range(c(cov(Y),1/(r+hi_wt-1)*sum_YYt)))
#   # image.plot(1/(r+hi_wt-1)*sum_YYt, zlim = range(c(cov(Y),1/(r+hi_wt-1)*sum_YYt)))  
#   Sigma_hat <- 1/(r+hi_wt-1) * sum_YYt
# } else {
#   Sigma_hat <- cov(Y)
# }
# 
# if(taper_cov) Sigma_hat <- Sigma_hat * bohman(plgp:::distance(x), tau=tau_b)
# 
# par(mfrow=c(1,2))
# matplot(x,t(Y),type="l",col="gray")
# lines(x,y_hi,col="blue",lwd=2)
# lines(x,y_avg,col="red",lwd=2)
# legend("bottomright", legend = c("hi","avg"), lty=1, col=c("blue","red"))
# image.plot(Sigma_hat/nrunn, main = "input as sigma_hat")

####################
# run for sim data #
####################

fitcov <- fit_two_layer_SW(x = x, y = y, nmcmc = nmcmc, cov = cov_fn, vecchia = vecchia, pmx = pmx, Sigma_hat = Lam_zi)
                            #settings = list(alpha =list(theta_w=1000), beta=list(theta_w=.0001/1000)))
# plot(fitcov)
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
plot.true.combo(fitcov)
lines(x, ((log10(y_hi)-temp_lm)-mean(y_un))/sd(y_un), col="red")
lines(x, ((log10(y_pt)-temp_lm)-mean(y_un))/sd(y_un), col="orange")
lines(x, ((log10(mu_z)-temp_lm)-mean(y_un))/sd(y_un), col="blue",lty=2)
lines(x, ((log10(y_lo_avg)-temp_lm)-mean(y_un))/sd(y_un))

if(!taper_cov) plot.true.combo(fitcov, S_e = fitcov$Sigma_hat*bohman(plgp:::distance(x),tau = tau_b))
lines(x, ((log10(y_hi)-temp_lm)-mean(y_un))/sd(y_un), col="red")
lines(x, ((log10(y_pt)-temp_lm)-mean(y_un))/sd(y_un), col="orange")
lines(x, ((log10(mu_z)-temp_lm)-mean(y_un))/sd(y_un), col="blue",lty=2)
lines(x, ((log10(y_lo_avg)-temp_lm)-mean(y_un))/sd(y_un))

if(!one_layer) plot.warp(fitcov, ref.scale = 1)

dev.off()

save.image(file = paste0("rda/1D_combo_",nmcmc,"_",one_layer,if(pmx){"_pmx"},if(vecchia){"_vec"},".rda"))
