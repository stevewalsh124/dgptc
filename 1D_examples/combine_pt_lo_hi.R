# Combine Perturbation Theory, Low Res, and Hi Res runs
# See "combine_pert_lo_hi.pdf" for modeling details
# real data study

vecchia <- F
pmx <- F
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
library(zoo)

# # No need to taper here! The covariance matrices are already diagonal.
# # Taper the covariance matrix before the MCMC fit?
# taper_cov <- F

# Smooth the precision info so it's not a step function for each data type (pt, lo, hi)?
smooth_precs <- T
if(smooth_precs) k_sm <- 10 #rolling mean uses k numbers

# Adjust effective variances for diagonal matrix
adj_var <- 1^2

# Do a kriging step?
krig <- F

ncores <- 2
tolpower <- -10

cov_fn <- "matern"#"exp2"#

tau_b <- 1
nrun <- 16 # number of low res runs; to adjust precision for average
nmcmc <- 6000
nburn <- 1000
kth <- 5

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

pdf(paste0("pdf/compare_pt_wts",nmcmc,"_",one_layer,adj_var,if(smooth_precs){paste0("_sp",k_sm)},
           if(pmx){"_pmx"},if(vecchia){"_vec"},".pdf"))
# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")

n <- length(k)

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
     ylim = range(c(log10(y_pt[index_list$pert.ix]), 
                    log10(y_lo_avg[index_list$lowres.ix]), 
                    log10(y_hi[index_list$highres.ix]))))
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

# subtract an overall average based on all 111*16 runs
x2 <- x^2
# temp_lm <- lm(log10(mu_z) ~ x + x2)$fitted.values
load(paste0("rda/avg_of_wt_avg",if(smooth_precs){paste0("_",k_sm)},".rda"))
temp_lm <- avg_loess
y_un <- c(log10(mu_z) - temp_lm)
y <- (y_un - mean(y_un))/sd(y_un)
lines(x, temp_lm, col="green")
legend("bottomleft",legend=c("pert","hi","lo","wt avg","lm"), 
       col = c("orange","red","black","blue","green"), lty=c(1,1,1,2,1))
lines(x[index_list$pert.ix], rep(2.35, length(index_list$pert.ix)))
lines(x[index_list$lowres.ix], rep(2.5, length(index_list$lowres.ix)), lty=2)
lines(x[index_list$highres.ix], rep(2.65, length(index_list$highres.ix)), lty=3)

# recompute upon adjusting for scaling (by var(y_un)) of log10(P(k))
Lam_lo <- diag(precs_lo)*var(y_un)/adj_var
Lam_hi <- diag(precs_hi)*var(y_un)/adj_var
Lam_pt <- diag(precs_pt)*var(y_un)
Lam_z <- Lam_pt+Lam_lo+Lam_hi
Lam_zi <- solve(Lam_z)

####################
# run for sim data #
####################
if(one_layer){
  fitcov <- fit_two_layer_SW(x = x, y = y, nmcmc = nmcmc, cov = cov_fn, 
                             vecchia = vecchia, Sigma_hat = Lam_zi)
  #settings = list(alpha =list(theta_w=1000), beta=list(theta_w=.0001/1000)))
} else {
  fitcov <- fit_two_layer_SW(x = x, y = y, nmcmc = nmcmc, cov = cov_fn, 
                             vecchia = vecchia, pmx = pmx, Sigma_hat = Lam_zi)
  #settings = list(alpha =list(theta_w=1000), beta=list(theta_w=.0001/1000)))
}

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
fitcov <- est.true.combo(fitcov)

plot.true.combo(fitcov)
lines(x, ((log10(y_hi)-temp_lm)-mean(y_un))/sd(y_un), col="red")
lines(x, ((log10(y_pt)-temp_lm)-mean(y_un))/sd(y_un), col="orange")
lines(x, ((log10(mu_z)-temp_lm)-mean(y_un))/sd(y_un), col="blue",lty=2)
lines(x, ((log10(y_lo_avg)-temp_lm)-mean(y_un))/sd(y_un))

plot.true.combo(fitcov, xlim = range(fitcov$x[index_list$lowres.ix]), 
                ylim=range(fitcov$lbb[index_list$lowres.ix],fitcov$ubb[index_list$lowres.ix]),
                legend.loc = "bottomright")
lines(x, ((log10(y_hi)-temp_lm)-mean(y_un))/sd(y_un), col="red")
lines(x, ((log10(y_pt)-temp_lm)-mean(y_un))/sd(y_un), col="orange")
lines(x, ((log10(mu_z)-temp_lm)-mean(y_un))/sd(y_un), col="blue",lty=2)
lines(x, ((log10(y_lo_avg)-temp_lm)-mean(y_un))/sd(y_un))

if(!one_layer) plot.warp(fitcov, ref.scale = 1)

dev.off()

save.image(file = paste0("rda/1D_combo_",nmcmc,"_",one_layer,adj_var,if(smooth_precs){paste0("_sp",k_sm)},
                         if(pmx){"_pmx"},if(vecchia){"_vec"},".rda"))
