# real data study

one_layer <- T
if(one_layer) {source("logl_cov_1L.R")} else {source("logl_cov.R")}

source("matrix.Moore.Penrose.R")
source("plot_fns.R") #plot.krig, plot.true, plot.warp

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
precs <- prec_lowres[index_list$lowres.ix]

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

bohman <- function(t, tau = 0.25){
  boh <- (1-t/tau)*cos(pi*t/tau)+sin(pi*t/tau)/pi
  boh <- ifelse(t>=tau, 0, boh)
}

cov_fn <- "matern"#"exp2"#

tau_b <- 1
nrun <- 16
nmcmc <- 7500
nburn <- 1500
kth <- 2

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

pdf(paste0("pdf/realstudydgpact_",nmcmc,"_",nrun,"_",
           if(taper_cov){paste0("tpr")},tau_b,
           if(use_hi){"_uh"},
           if(one_layer){"_1L"},"_",cov_fn,tolpower,".pdf"))

step <- 499

# Model 1, choose from 000-111
i <- 1

# first col is k, 2nd is linear pert theory, 3:18 is lowres, 19 is hires
pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# subtract a linear model
temp_avg <- rowMeans(log10(pk2[index_list$lowres.ix,3:18]))
temp_lm <- lm(temp_avg ~ log10(pk2[index_list$lowres.ix,1]))$fitted.values
Y <- t(apply(t(log10(pk2[index_list$lowres.ix,3:18])),1,function(x){x-temp_lm}))

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
y <- log10(pk2[index_list$lowres.ix, bte]) - temp_lm
y_lra <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]) - temp_lm)
y_hi <- log10(pk2[index_list$lowres.ix, 19]) - temp_lm

# hi res precs on the low res index (for comparison)
precs_hi <- prec_highres[index_list$lowres.ix]
if(use_hi){
  # Weighted average of low res runs and 
  hi_wt <- unique(prec_highres/prec_lowres)[1]
  y_avg <- (hi_wt*y_hi+nrun*y_lra)/(hi_wt+nrun)
  nrunn <- nrun + hi_wt
} else {
  y_avg <- y_lra
  nrunn <- nrun
}

# wavenumber is X, a particular lowres run in Y
mean_sz <- mean(y_avg)
sd_sz <- sd(y_avg)
x <- log10(k[index_list$lowres.ix])
x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean_sz)/sd_sz
y_avg <- (y_avg - mean_sz)/sd_sz
y_lra <- (y_lra - mean_sz)/sd_sz
y_hi <- (y_hi - mean_sz)/sd_sz
for (i in 1:nrow(Y)) { Y[i,] <- (Y[i,] - mean_sz)/sd_sz  }

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

if(taper_cov) Sigma_hat <- Sigma_hat * bohman(plgp:::distance(x), tau=tau_b)

par(mfrow=c(1,2))
matplot(x,t(Y),type="l",col="gray")
lines(x,y_hi,col="blue",lwd=2)
lines(x,y_avg,col="red",lwd=2)
legend("bottomright", legend = c("hi","avg"), lty=1, col=c("blue","red"))
image.plot(Sigma_hat/nrunn, main = "input as sigma_hat")

####################
# run for sim data #
####################

fitcov <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = nmcmc, Sigma_hat = Sigma_hat/nrunn, cov = cov_fn)
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

par(mfrow=c(1,2))
plot.true(fitcov)
if(!taper_cov) plot.true(fitcov, S_e = fitcov$Sigma_hat*bohman(plgp:::distance(x),tau = tau_b))

if(!one_layer) plot.warp(fitcov)

dev.off()