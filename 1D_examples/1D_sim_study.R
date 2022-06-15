# sim study

one_layer <- F
if(one_layer) {source("logl_cov_1L.R")} else {source("logl_cov.R")}

source("matrix.Moore.Penrose.R")
source("plot_fns.R")

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
precs <- prec_lowres[index_list$lowres.ix]

library(MASS) #ginv
library(fields) #image.plot
library(mvtnorm) #rmvnorm
library(plgp) #distance (which is squared distances)

# Should the true error's cov mtx be diagonal or dense?
true_diag <- F
# Use the true covariance matrix for sigma hat?
use_true_cov <- F
# Taper the covariance matrix?
taper_cov <- F

seed <- 1
for (seed in 1:10) {

cov_fn <- "matern"#"exp2"#

nrun <- 16
nmcmc <- 22500
nburn <- 2500
kth <- 4

pdf(paste0("pdf/simstudydgpact_",nmcmc,"_",nrun,
           if(true_diag){"_TD"},
           if(use_true_cov){"_UTC"},
           if(taper_cov){"_tpr"},
           if(one_layer){"_1L"},"_",cov_fn,"_",seed,".pdf"))

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

step <- 499
# Model 1, choose from 000-111
i <- 1

pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# hi res precs on the low res index (for comparison)
precs_hi <- prec_highres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
# y <- log10(pk2[index_list$lowres.ix, bte])
# Y <- log10(pk2[index_list$lowres.ix, 3:18])
x <- (x - min(x))/(max(x)-min(x))
# y <- (y - mean(y))/sd(y)
# y_avg <- rowMeans(Y)
# y_avg <- (y_avg - mean(y_avg))/sd(y_avg)

# computes squared distances
D <- plgp:::distance(x)

# lengthscale = 1
# marginal variance = 1
Sigma <- exp(-D) + diag(eps, length(x))

# first layer
set.seed(seed)
w <- x^3#c(rmvnorm(1, mean=x, sigma=Sigma))

# second layer
Dw <- plgp:::distance(c(w))
Sigma_s <- exp(-Dw/0.05) + diag(eps, length(w))
ytrue <- rmvnorm(1, sigma=Sigma_s) #mean=W,

#Plot each combination of layers/output
par(mfrow=c(1,3))
plot(x, w, type="l"); abline(0, 1, col=2, lty=2)
plot(w, ytrue, type = "l"); abline(h=0, col=2, lty=2)
plot(x, ytrue, type="l"); abline(h=0, col=2, lty=2)

# lmf <- lm(y_avg ~ x)$fitted.values
# ytrue <- y_avg - lmf
# ytrue <- (ytrue - mean(ytrue))/sd(ytrue)
# plot(ytrue, type="l", main = "true y")

# generate 16 low res runs with iid noise (to start)
# Gaussian cov fn with theta=1, tau2=1
if(true_diag){
  Cov_true <- diag(5000/precs^0.8)
} else {
  Sig_mat <- exp(-plgp:::distance(x/.05)) + diag(sqrt(.Machine$double.eps),length(x)) 
  A <- diag(2500/precs^0.8)
  Cov_true <- A %*% Sig_mat %*% A
}

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
  Sigma_hat <- cov(Y_sim)
}

if(taper_cov){
  Sigma_taper <- matrix(0, nrow(Sigma_hat), ncol(Sigma_hat))
  for (i in 1:nrow(Sigma_hat)) {
    for (j in (i-5):(i+5)) {
      if(j <= 0 || j > nrow(Sigma_hat)) next
      Sigma_taper[i,j] <- Sigma_hat[i,j]
    }
  }
  Sigma_hat <- Sigma_taper
}

par(mfrow=c(1,2))
image.plot(Cov_true, main = "Cov_true", zlim = range(c(Cov_true,Sigma_hat)))
image.plot(Sigma_hat, main = "input as sigma_hat", zlim = range(c(Cov_true,Sigma_hat)))

####################
# run for sim data #
####################

fitcov <- fit_two_layer_SW(x = x, y = colMeans(Y_sim), nmcmc = nmcmc, 
                           Sigma_hat = Sigma_hat/nrun, cov = cov_fn)
# plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
plot(fitcov)
fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=2, precs_pred = 1/diag(Cov_true_pred))

zz <- fitcov$mean#-fitcov$mean

par(mfrow=c(1,1))
plot(fitcov$x_new, zz, type="l", col="blue", lwd=1.5,
     ylim = range(c(#zz-2*sqrt(fitcov$s2_smooth*mean(fitcov$tau2)*mean(fitcov$g)),
                    #zz+2*sqrt(fitcov$s2_smooth*mean(fitcov$tau2)*mean(fitcov$g)),
                    zz-2*sqrt(fitcov$s2*mean(fitcov$tau2)),
                    zz+2*sqrt(fitcov$s2*mean(fitcov$tau2)),
                    zz-2*sqrt(1/precs_pred),
                    zz+2*sqrt(1/precs_pred), Y_sim)), main = "fitcov")
for (i in 1:nrow(Y_sim))  lines(c(x), Y_sim[i,], col="gray")#-y_avg)
lines(fitcov$x_new, zz, col="blue", lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov$x_new, zz-2*sqrt(1/(nrun*precs_pred)), col=3, lwd=1.5)
lines(fitcov$x_new, zz+2*sqrt(1/(nrun*precs_pred)), col=3, lwd=1.5)
lines(fitcov$x, ytrue, lwd=1.5, lty=2)
legend("topright", col=c("blue",2:3,1), legend = c("mean","UQ", "precs (low res)","true S"), 
       lty=c(1,1,1,2), lwd=1.5)
# save(fitcov, file = paste0("rda/sim/fitcov",nmcmc,if(one_layer){"_1L"},"_lm.rda"))

v <- fitcov$v

# Sigma_epsilon; already adjusted by 16 and sigma2_yavg when rescaling data
# generalized (Moore-Penrose) inverse
# for (ep in c(0,1e-10,1e-5,1e-1)) {

S_ei <- matrix.Moore.Penrose(fitcov$Sigma_hat)# + diag(rep(1e-7,161)))
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
  if(v==999){
    S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
  } else {
    S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
  }
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

plot(fitcov$x, fitcov$y, type="n",
     ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), main = paste0("fitcov, ",nrun,"-sample cov mtx"))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

# Using law of total expectation and variance instead
mu_bar <- rowMeans(Ms)
C_bar <- matrix(rowMeans(Cs),  length(fitcov$x),  length(fitcov$x))
cov_mut <- cov(t(Ms))

par(mfrow=c(1,2))
plot(x, mu_bar, type="n", main = "law of total E & V; MPI version", col="blue",
     ylim=range(Y_sim))
for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(x, mu_bar, col="blue")
lines(x, ytrue, col="red")
lines(x, fitcov$y, lty=2)
lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)), col="blue")
lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)), col="blue")
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

# true error covariance, estimated prior cov
S_ei <- solve(Cov_true/16)# + diag(rep(1e-7,161)))
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
  if(v==999){
    S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
  } else {
    S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
  }
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

plot(fitcov$x, fitcov$y, type="n",
     ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), main = paste0("est prior matrix, true cov mtx"))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

# est error cov (tapered), true prior cov
sighat_temp <- cov(Y_sim)*abs(cor(Y_sim)>0.1)
sighattemp_inv <- matrix.Moore.Penrose2(sighat_temp/16)
Ctrue <- solve(solve(Sigma_s) + sighattemp_inv)
Ctrue <- (Ctrue+t(Ctrue))/2
Mtrue <- Ctrue %*% sighattemp_inv %*% fitcov$y
true_samp <- matrix(NA, nrow(Ss), 1000)
for(i in 1:1000) true_samp[,i] <- Mtrue + matrix.sqrt(Ctrue)%*%rnorm(length(fitcov$y))

m <- rowMeans(true_samp) #- y_avg
ub <- apply(true_samp, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(true_samp, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(true_samp, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(true_samp, 1, function(x){quantile(x,0.005)}) #- y_avg

plot(fitcov$x, fitcov$y, type="n",ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), 
     main = paste0(nrun,"-sample cov mtx tapered,\n but true prior"))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

# est error cov (NO taper), true prior cov
sighat_temp <- Sigma_hat
sighattemp_inv <- matrix.Moore.Penrose2(sighat_temp/16)
Ctrue <- solve(solve(Sigma_s) + sighattemp_inv)
Ctrue <- (Ctrue+t(Ctrue))/2
Mtrue <- Ctrue %*% sighattemp_inv %*% fitcov$y
true_samp <- matrix(NA, nrow(Ss), 1000)
for(i in 1:1000) true_samp[,i] <- Mtrue + matrix.sqrt(Ctrue)%*%rnorm(length(fitcov$y))

m <- rowMeans(true_samp) #- y_avg
ub <- apply(true_samp, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(true_samp, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(true_samp, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(true_samp, 1, function(x){quantile(x,0.005)}) #- y_avg

plot(fitcov$x, fitcov$y, type="n",ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), 
     main = paste0(nrun,"-sample cov mtx NO taper,\n but true prior"))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

# generate true uncertainty
Ctrue <- solve(solve(Sigma_s) + solve(Cov_true/16))
Ctrue <- (Ctrue+t(Ctrue))/2
Mtrue <- Ctrue %*% solve(Cov_true/16) %*% fitcov$y
true_samp <- matrix(NA, nrow(Ss), 2500)
for(i in 1:2500) true_samp[,i] <- Mtrue + matrix.sqrt(Ctrue)%*%rnorm(length(fitcov$y))

m <- rowMeans(true_samp) #- y_avg
ub <- apply(true_samp, 1, function(x){quantile(x,0.975)}) #- y_avg
lb <- apply(true_samp, 1, function(x){quantile(x,0.025)}) #- y_avg
ubb <- apply(true_samp, 1, function(x){quantile(x,0.995)}) #- y_avg
lbb <- apply(true_samp, 1, function(x){quantile(x,0.005)}) #- y_avg

emp_cover <- round(mean(ytrue > lb & ytrue < ub),3)
emp_cover99 <- round(mean(ytrue > lbb & ytrue < ubb),3)

plot(fitcov$x, fitcov$y, type="n",ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), 
     main = paste0("fitcov, ",nrun,"-samples, \nboth mtxs true\n",
                   emp_cover, " ", emp_cover99))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

bohman <- function(t, tau = 0.25){
  boh <- (1-t/tau)*cos(pi*t/tau)+sin(pi*t/tau)/pi
  boh <- ifelse(t>=tau, 0, boh)
}

S_e <- Sigma_hat/nrun * bohman(plgp:::distance(x), tau=0.1)
S_ei <- matrix.Moore.Penrose2(S_e)
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
  if(v==999){
    S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
  } else {
    S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
  }
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

emp_cover <- round(mean(ytrue > lb & ytrue < ub),3)
emp_cover99 <- round(mean(ytrue > lbb & ytrue < ubb),3)

plot(fitcov$x, fitcov$y, type="n",
     ylim = range(c(m, lb, ub, lbb, ubb, Y_sim)), main = paste0("est both, taper error w bohman\n",
                                                                emp_cover, " ", emp_cover99))

for (i in 1:nrun) lines(fitcov$x, Y_sim[i,], col="gray")
lines(fitcov$x, ytrue, lwd=1.5, col="red")
lines(fitcov$x, fitcov$y, lwd=1.5, lty=2)
lines(fitcov$x, m , col="blue")
lines(fitcov$x, lb , col="blue", lty=2)
lines(fitcov$x, ub , col="blue", lty=2)
lines(fitcov$x, lbb , col="darkblue", lty=2)
lines(fitcov$x, ubb , col="darkblue", lty=2)
legend("bottomright", legend = c("true", "sample avg", "UQ"), 
       col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))

if(!one_layer) plot.warp(fitcov)

dev.off()

}
