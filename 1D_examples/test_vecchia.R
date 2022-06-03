one_layer <- F
if(one_layer) {source("logl_cov_1L.R")} else {source("logl_cov.R")}

source("matrix.Moore.Penrose.R")

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
precs <- prec_lowres[index_list$lowres.ix]

library(MASS) #ginv
library(fields) #image.plot

do_diag <- T

nmcmc <- 13101#3
nburn <- 1101#3
kth <- 1

pdf(paste0("pdf/uq_sim_truth_fixtau2_mean_lmybar_",nmcmc,if(one_layer){"_1L"},".pdf"))

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

step <- 499
# Model 1, choose from 000-111
i <- 1

# gather all cov mtx for the 111 models
covYms <- covYts <- list()
YYm <- YYt <- YY <- c()
mses <- ppk2 <- ppk2a <- ppk2b <- c()

for (i in 1:111) {
  # first col is k, 2nd is linear pert theory, 3:18 is lowres, 19 is hires
  pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                           if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
  ppk2 <- rbind(ppk2,t(log10(pk2[index_list$lowres.ix,3:18])))
  
  temp_avg <- rowMeans(log10(pk2[index_list$lowres.ix,3:18]))
  temp_pk2a <- t(apply(t(log10(pk2[index_list$lowres.ix,3:18])),1,function(x){x-temp_avg}))
  ppk2a <- rbind(ppk2a,temp_pk2a)
  
  temp_lm <- lm(temp_avg ~ log10(pk2[index_list$lowres.ix,1]))$fitted.values
  temp_pk2a <- t(apply(t(log10(pk2[index_list$lowres.ix,3:18])),1,function(x){x-temp_lm}))
  ppk2b <- rbind(ppk2b,temp_pk2a)
}

for(i in 1:111) mses[i] <- mean((cov(ppk2[1:16,]) -  cov(ppk2[16*(i-1)+(1:16),]))^2)
mses_to_use <- order(mses)[1:20]

for (i in length(mses_to_use):1) {
  ind <- mses_to_use[i]
  pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                           if(ind<100){"0"},if(ind<10){"0"},ind,"_test.dat"))
  
  # hi res precs on the low res index (for comparison)
  precs_hi <- prec_highres[index_list$lowres.ix]
  
  # wavenumber is X, a particular lowres run in Y
  x <- log10(k[index_list$lowres.ix])
  y <- log10(pk2[index_list$lowres.ix, bte])
  y_avg <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]))
  y_hi <- log10(pk2[index_list$lowres.ix, 19])
  
  # get sample covariance matrix for the 16 low res runs
  Y <- t(log10(pk2[index_list$lowres.ix, 3:18]))
  YY <- rbind(YY,Y)
  Yt <- t(apply(Y,1,function(x){x-y_avg}))
  covYts[[ind]] <- cov(Yt)
  YYt <- rbind(YYt, Yt)
  
  ybarlm <- lm(y_avg ~ x)$fitted.values
  Ym <- t(apply(Y,1,function(x){x-ybarlm}))
  covYms[[ind]] <- cov(Ym)
  YYm <- rbind(YYm, Ym)
  
}

# Plot with Y bar subtracted
plot(x,Y[1,]-y_avg, type="l", main = paste("Y-y_avg",i), ylim = range(Yt), col="gray")
for (k in 2:16) lines(x,Y[k,]-y_avg,col="gray")
lines(x,y_avg-y_avg, col="red")
legend("topright", legend=c("avg"), col=c("red"), lty=1)
Ym <- t(apply(Y,1,function(x){x-ybarlm}))
ym_avg <- colMeans(Ym)

# subtract a linear fit for a better prior mean
plot(x,Ym[1,], type="l", main = "Y minus lm for avg", ylim = range(Ym), col="gray")
for (k in 2:16) lines(x,Ym[k,],col="gray")
lines(x, ybarlm - ybarlm, col="green")
lines(x, ym_avg, col="red", lwd=2)
legend("bottomright", legend=c("avg","lm"), col=c("red","green"), lty=1)

# save var(y_avg) to adjust prec values in the final fit call
sigma2_yavgm <- sd(Ym)^2

x <- (x - min(x))/(max(x)-min(x))
# y <- (y - mean(y))/sd(y)
# Y <- (Y - mean(Y))/sd(Y)
Ym <- (Ym - mean(Ym))/sd(Ym)
ym_avg <- colMeans(Ym)
y_hi <- (y_hi - mean(y_hi))/sd(y_hi)

# get average cov mtx
covYms_comb <- Reduce("+", covYms)/length(covYms)
covYts_comb <- Reduce("+", covYts)/length(covYts)
cov_YYm <- cov(YYm)
cov_YYt <- cov(YYt)

# obtain sample cov mat after rescaling
covYm <- cov(Ym)
# covYY <- cov(YY)
# plot(log10(diag(covYY)))

# get predictions for xx and have corresponding prec info
xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(precs) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))

#######################
# run for actual data #
#######################

Yadj <- Ym#apply(Y,1,function(x){x-y_avg})

if(do_diag){
  fitdiag <- fit_two_layer_SW(x = x, y = ym_avg, nmcmc = nmcmc, Sigma_hat = diag(1/(precs*16*sigma2_yavgm)))
  # plot(fitdiag)
  fitdiag <- trim_SW(fitdiag, nburn, kth)
  plot(fitdiag)
  fitdiag <- predict.dgp2_SW(object = fitdiag, xx, cores=2, precs_pred = precs_pred)
  
  zz <- fitdiag$mean#-fitdiag$mean
  
  par(mfrow=c(1,1))
  plot(fitdiag$x_new, zz, type="l", col="blue", lwd=1.5,
       ylim = range(c(zz-2*sqrt(fitdiag$s2_smooth*mean(fitdiag$tau2)*mean(fitdiag$g)),
                      zz+2*sqrt(fitdiag$s2_smooth*mean(fitdiag$tau2)*mean(fitdiag$g)),
                      zz-2*sqrt(fitdiag$s2*mean(fitdiag$tau2)),
                      zz+2*sqrt(fitdiag$s2*mean(fitdiag$tau2)),
                      zz-2*sqrt(1/precs_pred),
                      zz+2*sqrt(1/precs_pred), Yadj)), main = "fitdiag")
  for (i in 1:nrow(Y))  lines(c(x), Ym[i,])#-y_avg)
  lines(fitdiag$x_new, zz, col="blue", lwd=1.5)
  lines(fitdiag$x_new, zz-2*sqrt(fitdiag$s2_smooth), col=2, lwd=1.5)
  lines(fitdiag$x_new, zz+2*sqrt(fitdiag$s2_smooth), col=2, lwd=1.5)
  lines(fitdiag$x_new, zz-2*sqrt(fitdiag$s2), col=2, lwd=1.5)
  lines(fitdiag$x_new, zz+2*sqrt(fitdiag$s2), col=2, lwd=1.5)
  lines(fitdiag$x_new, zz-2*sqrt(1/precs_pred), col=3, lwd=1.5)
  lines(fitdiag$x_new, zz+2*sqrt(1/precs_pred), col=3, lwd=1.5)
  lines(fitdiag$x_new, zz-2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
  lines(fitdiag$x_new, zz+2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
  legend("topright", col=c("blue",2:3), legend = c("mean","UQ", "precs (low res)"), lty=1, lwd=1.5)
  save(fitdiag, file = paste0("rda/corr_err/fitdiag",nmcmc,if(one_layer){"_1L"},"_lm.rda"))
}

xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(1/diag(covYm)) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))


fitcov <- fit_two_layer_SW(x = x, y = ym_avg, nmcmc = nmcmc, Sigma_hat = cov(Yt)/16)#, cov = "exp2")
# plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
plot(fitcov)
fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=2, precs_pred = precs_pred)

zz <- fitcov$mean

par(mfrow=c(1,1))
plot(fitcov$x_new, zz, type="l", col="blue", lwd=1.5, 
     ylim = range(c(zz-2*sqrt(fitcov$s2_smooth),
                    zz+2*sqrt(fitcov$s2_smooth),
                    zz-2*sqrt(fitcov$s2),
                    zz+2*sqrt(fitcov$s2),
                    zz-2*sqrt(1/precs_pred),
                    zz+2*sqrt(1/precs_pred), Yadj)), main = "fitcov")
for (i in 1:nrow(Y))  lines(c(x), Ym[i,])#-y_avg)
lines(fitcov$x_new, zz, col="blue", lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
legend("topright", col=c("blue",2:3), legend = c("mean","UQ", "precs (low res)"), lty=1, lwd=1.5)
save(fitcov, file = paste0("rda/corr_err/fitcov",nmcmc,if(one_layer){"_1L"},"_lm.rda"))


fitcov_fr <- fit_two_layer_SW(x = x, y = ym_avg, nmcmc = nmcmc, Sigma_hat = cov_YYt/16)#, cov = "exp2")
# plot(fitcov_fr)
fitcov_fr <- trim_SW(fitcov_fr, nburn, kth)
plot(fitcov_fr)
fitcov_fr <- predict.dgp2_SW(object = fitcov_fr, xx, cores=2, precs_pred = precs_pred)

zz <- fitcov$mean#-fitcov$mean

par(mfrow=c(1,1))
plot(fitcov_fr$x_new, zz, type="l", col="blue", lwd=1.5, 
     ylim = range(c(zz-2*sqrt(fitcov_fr$s2_smooth),
                    zz+2*sqrt(fitcov_fr$s2_smooth),
                    zz-2*sqrt(fitcov_fr$s2),
                    zz+2*sqrt(fitcov_fr$s2),
                    zz-2*sqrt(1/precs_pred),
                    zz+2*sqrt(1/precs_pred), Ym)), main = "fitcov_fr")
for (i in 1:nrow(Y))  lines(c(x), Ym[i,])#-y_avg)
lines(fitcov_fr$x_new, zz, col="blue", lwd=1.5)
lines(fitcov_fr$x_new, zz-2*sqrt(fitcov_fr$s2_smooth), col=2, lwd=1.5)
lines(fitcov_fr$x_new, zz+2*sqrt(fitcov_fr$s2_smooth), col=2, lwd=1.5)
lines(fitcov_fr$x_new, zz-2*sqrt(fitcov_fr$s2), col=2, lwd=1.5)
lines(fitcov_fr$x_new, zz+2*sqrt(fitcov_fr$s2), col=2, lwd=1.5)
lines(fitcov_fr$x_new, zz-2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov_fr$x_new, zz+2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov_fr$x_new, zz-2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
lines(fitcov_fr$x_new, zz+2*sqrt(1/(16*precs_pred)), col=3, lwd=1.5)
legend("topright", col=c("blue",2:3), legend = c("mean","UQ", "precs (low res)"), lty=1, lwd=1.5)
save(fitcov_fr, file = paste0("rda/corr_err/fitcov_fr",nmcmc,if(one_layer){"_1L"},"_lm.rda"))


# UQ for truth S by Bayes Theorem:
# S | \bar{Y} \sim(m, C) where
# C^{-1} = (Sigma_s(x,x'))^{-1} + 16*Sigma_epsilon^{-}
# m = C * 16*Sigma_epsilon^{-} \bar{Y}
# In particular, using the MCMC output for the unknown parameters in Sigma_s(x,x'), 
# we can use the distribution [S | \bar{Y}] above to simulate a sample from posterior of S. 

v <- fitcov$v

# Sigma_epsilon; already adjusted by 16 and sigma2_yavg when rescaling data
# generalized (Moore-Penrose) inverse
# for (ep in c(0,1e-10,1e-5,1e-1)) {

S_ei <- matrix.Moore.Penrose(fitcov$Sigma_hat)# + diag(rep(ep,161)))
S_ei <- (S_ei+t(S_ei))/2

Cs <- matrix(NA, length(fitcov$x)^2, fitcov$nmcmc)
Ss <- St <- Sw <- Sx <- Ms <- matrix(NA, length(fitcov$x), fitcov$nmcmc)
for (i in 1:fitcov$nmcmc){
  if( i %% 2500 == 0) print(i)
  theta <- fitcov$theta_y[i]
  # tau2 <- fitcov$tau2[i]
  # t2S_ei <- 1/tau2 * S_ei
  g <- fitcov$g[i]
  W <- fitcov$w[[i]]
  dw <- sq_dist(W)
  
  S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
  S_si <- (S_si+t(S_si))/2
  C <- Cs[,i] <- solve(S_si + S_ei)
  C <- (C+t(C))/2
  M <- Ms[,i] <- C %*% S_ei %*% fitcov$y
  if(i %% 1000 ==0) print(range(M))
  Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fitcov$y))
  # St[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "eigen")
  Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
  # Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
}

m <- rowMeans(Ss) #- y_avg
ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg

plot(fitcov$x, fitcov$y, type="l",
     ylim = range(c(m, lb, ub, lbb, ubb, Yadj, Ym)), main = paste0("fitcov, 16-sample cov mtx"))

for (i in 1:16) lines(fitcov$x, Ym[i,], col="gray")
lines(fitcov$x, fitcov$y - fitcov$y, lwd=1.5)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)

mu_bar <- rowMeans(Ms)
C_bar <- matrix(rowMeans(Cs),  length(fitcov$x),  length(fitcov$x))
cov_mut <- cov(t(Ms))

# Using law of total expectation and variance instead
plot(x, mu_bar, type="l", main = "law of total E & V; MPI version")
lines(x, ym_avg, col="blue")
lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)))
lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)))
# }

# invertible cov mtx from 111*16 low res runs
S_ei <- solve(fitcov_fr$Sigma_hat)
S_ei <- (S_ei+t(S_ei))/2

Cs <- matrix(NA, length(fitcov$x)^2, fitcov$nmcmc)
Ss <- St <- Sw <- Sx <- Ms <- matrix(NA, length(fitcov$x), fitcov$nmcmc)
for (i in 10000:fitcov_fr$nmcmc){
  print(i)#if( i %% 2500 == 0) print(i)
  theta <- fitcov_fr$theta_y[i]
  tau2 <- fitcov_fr$tau2[i]
  t2S_ei <- 1/tau2 * S_ei
  g <- fitcov_fr$g[i]
  W <- fitcov_fr$w[[i]]
  dw <- sq_dist(W)
  
  S_si <- 1/tau2 * solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v)))
  S_si <- (S_si+t(S_si))/2
  C <- solve(S_si + t2S_ei)
  C <- Cs[,i] <- (C + t(C))/2
  M <- Ms[,i] <- C %*% t2S_ei %*% fitcov_fr$y
  if(i %% 1000 ==0) print(range(M))
  Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fitcov_fr$y))
  St[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "eigen")
  Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
  Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
}

m <- rowMeans(Ss) #- y_avg
ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg

plot(fitcov_fr$x, fitcov_fr$y, type="l",
     ylim = range(c(m, lb, ub, Yadj)), main = "fitcov_fr, 1776-sample cov mtx")
for (i in 1:16) lines(fitcov_fr$x, Ym[i,], col="gray")
lines(fitcov_fr$x, fitcov_fr$y, lwd=1.5)
abline(h=0, col="green")
lines(fitcov_fr$x, m , col="blue")
lines(fitcov_fr$x, lb , col="blue", lty=2)
lines(fitcov_fr$x, ub , col="blue", lty=2)
lines(fitcov_fr$x, lbb , col="darkblue", lty=2)
lines(fitcov_fr$x, ubb , col="darkblue", lty=2)

mu_bar <- rowMeans(Ms)
C_bar <- matrix(rowMeans(Cs),  length(fitcov$x),  length(fitcov$x))
cov_mut <- cov(t(Ms))

# Using law of total expectation and variance instead
plot(x, mu_bar, type="l", main = "law of total E & V; full rank version")
lines(x, ym_avg, col="blue")
lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)))
lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)))

# invertible cov mtx from 111*16 low res runs
S_ei <- solve(cov(ppk2a)/16)
S_ei <- (S_ei+t(S_ei))/2
Cs <- matrix(NA, length(fitcov$x)^2, fitcov$nmcmc)
Ss <- St <- Sw <- Sx <- Ms <- matrix(NA, length(fitcov$x), fitcov$nmcmc)
for (i in 1:fitcov_fr$nmcmc){
  if( i %% 2500 == 0) print(i)
  theta <- fitcov_fr$theta_y[i]
  tau2 <- fitcov_fr$tau2[i]
  t2S_ei <- 1/tau2 * S_ei
  g <- fitcov_fr$g[i]
  W <- fitcov_fr$w[[i]]
  dw <- sq_dist(W)
  
  S_si <- 1/tau2 * solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v)))
  S_si <- (S_si+t(S_si))/2
  C <- solve(S_si + t2S_ei)
  C <- Cs[,i] <- (C + t(C))/2
  M <- Ms[,i] <- C %*% t2S_ei %*% fitcov_fr$y
  if(i %% 1000 ==0) print(range(M))
  Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fitcov_fr$y))
  St[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "eigen")
  Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
  Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
}

mu_bar <- rowMeans(Ms)
C_bar <- matrix(rowMeans(Cs),  length(fitcov$x),  length(fitcov$x))
cov_mut <- cov(t(Ms))

# Using law of total expectation and variance instead
plot(x, mu_bar, type="l", main = "law of total E & V; ppk2 version")
lines(x, ym_avg, col="blue")
lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)))
lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)))

m <- rowMeans(Ss) #- y_avg
ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg

plot(fitcov_fr$x, fitcov_fr$y, type="l",
     ylim = range(c(m, lb, ub, Yadj)), main = "fitcov_fr, ppk2 version")
for (i in 1:16) lines(fitcov_fr$x, Ym[i,], col="gray")
lines(fitcov_fr$x, fitcov_fr$y, lwd=1.5)
abline(h=0, col="green")
lines(fitcov_fr$x, m , col="blue")
lines(fitcov_fr$x, lb , col="blue", lty=2)
lines(fitcov_fr$x, ub , col="blue", lty=2)
lines(fitcov_fr$x, lbb , col="darkblue", lty=2)
lines(fitcov_fr$x, ubb , col="darkblue", lty=2)

if(do_diag){
  # Plot the diagonal case
  S_ei <- solve(fitdiag$Sigma_hat)#diag(16*sigma2_yavgm*precs)
  
  Ss <- matrix(NA, length(fitdiag$x), fitdiag$nmcmc)
  for (i in 1:fitdiag$nmcmc){
    if( i %% 2500 == 0) print(i)
    theta <- fitdiag$theta_y[i]
    tau2 <- fitdiag$tau2[i]
    t2S_ei <- 1/tau2 * S_ei
    g <- fitdiag$g[i]
    W <- fitdiag$w[[i]]
    dw <- sq_dist(W)
    
    S_si <- 1/tau2 * solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v)))
    C <- solve(S_si + t2S_ei)
    M <- C %*% t2S_ei %*% fitdiag$y
    if(i %% 1000 ==0) print(range(M))
    Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fitdiag$y))
  }
  
  m <- rowMeans(Ss) #- y_avg
  ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
  lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
  ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
  lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg
  
  plot(fitdiag$x, fitdiag$y, type="l",
       ylim = range(c(Yadj, m, lb, ub)), main = "fitdiag, given prec mtx")
  
  for (i in 1:16) lines(fitdiag$x, Ym[i,], col="gray")
  lines(fitdiag$x, fitdiag$y - fitdiag$y, lwd=1.5)
  lines(fitdiag$x, m , col="blue")
  lines(fitdiag$x, lb , col="blue", lty=2)
  lines(fitdiag$x, ub , col="blue", lty=2)
  lines(fitdiag$x, lbb , col="darkblue", lty=2)
  lines(fitdiag$x, ubb , col="darkblue", lty=2)
  
  par(mfrow=c(2,2))
  image.plot(fitdiag$Sigma_hat, main = "fitdiag")
  image.plot(solve(fitdiag$Sigma_hat), main = "fitdiag inv")
  image.plot(matrix.Moore.Penrose(fitdiag$Sigma_hat), main = "fitdiag MPi")
  image.plot(matrix.Moore.Penrose(fitdiag$Sigma_hat)%*%fitdiag$Sigma_hat, main = "MPi * Sigma")
}

par(mfrow=c(2,2))
image.plot(fitcov$Sigma_hat, main = "fitcov")
image.plot(solve(fitcov$Sigma_hat), main = "fitcov inv")
image.plot(matrix.Moore.Penrose(fitcov$Sigma_hat), main = "fitcov MPi")
image.plot(matrix.Moore.Penrose(fitcov$Sigma_hat)%*%fitcov$Sigma_hat, main = "MPi * Sigma")


par(mfrow=c(2,2))
image.plot(fitcov_fr$Sigma_hat, main = "fitcov_fr")
image.plot(solve(fitcov_fr$Sigma_hat), main = "fitcov_fr inv")
image.plot(matrix.Moore.Penrose(fitcov_fr$Sigma_hat), main = "fitcov_fr MPi")
image.plot(matrix.Moore.Penrose(fitcov_fr$Sigma_hat)%*%fitcov_fr$Sigma_hat, main = "MPi * Sigma")

# # try another covariance matrix
# ppk2 <- c()
# for (i in 1:111) {
#   pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
#                            if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
#   mses <- ppk2 <- c()
#   for(i in 1:111) mses[i] <- mean((cov(ppk2[1:16,]) -  cov(ppk2[16*(i-1)+(1:16),]))^2)
#   
#   ppk2 <- rbind(ppk2,t(log10(pk2[index_list$lowres.ix,3:18])))
# }

par(mfrow=c(2,2))
image.plot(cov(ppk2[1:32,]), main = "ppk2 all batches")
image.plot(cov(ppk2[1:16,]), main = "ppk2 one batch")
image.plot(solve(cov(ppk2)), main = "ppk2 inv")
image.plot(matrix.Moore.Penrose(cov(ppk2)), main = "ppk2 MPi")

# pdf("pdf/cov_cor_eachbatch.pdf")
# par(mfrow=c(3,3))
# for(i in 1:111) image.plot(cov(ppk2[16*(i-1)+(1:16),]), main = paste("cov ppk2 batch",i))
# 
# par(mfrow=c(3,3))
# for(i in 1:111) image.plot(cor(ppk2[16*(i-1)+(1:16),]), main = paste("cor ppk2 batch",i))

dev.off()
