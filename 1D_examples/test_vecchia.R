source("logl_cov.R")

nmcmc <- 52500
nburn <- 2500
kth <- 5

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

step <- 499
for (i in 111:1) {
  
# i <- 1 # Model 1, choose from 000-111
pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
precs <- prec_lowres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
y <- log10(pk2[index_list$lowres.ix, bte])
y_avg <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]))
y_hi <- log10(pk2[index_list$lowres.ix, 19])

# get sample covariance matrix for the 16 low res runs
Y <- t(log10(pk2[index_list$lowres.ix, 3:18]))

# save var(y_avg) to adjust prec values in the final fit call
sigma2_y <- var(y)
sigma2_yavg <- var(y_avg)

x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean(y))/sd(y)
Y <- (Y - mean(Y))/sd(Y)
y_avg <- (y_avg - mean(y_avg))/sd(y_avg)
y_hi <- (y_hi - mean(y_hi))/sd(y_hi)

if(i == 111){ YY <- Y; plot(log10(diag(cov(Y))))} else {YY <- rbind(YY,Y)}
}
# obtain sample cov mat after rescaling
covY <- cov(Y)
covYY <- cov(YY)
plot(log10(diag(covYY)))

# get predictions for xx and have corresponding prec info
xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(precs) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))

#######################
# run for actual data #
#######################

fitdiag <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = nmcmc, Sigma_hat = diag(1/precs*16*sigma2_yavg))
# plot(fitdiag)
fitdiag <- trim_SW(fitdiag, nburn, kth)
plot(fitdiag)
fitdiag <- predict.dgp2_SW(object = fitdiag, xx, cores=6, precs_pred = precs_pred)

zz <- fitdiag$mean-fitdiag$mean
Yadj <- apply(Y,1,function(x){x-y_avg})

par(mfrow=c(1,1))
plot(fitdiag$x_new, zz, type="l", col="blue", lwd=1.5,
     ylim = range(c(zz-2*sqrt(fitdiag$s2_smooth*mean(fitdiag$tau2)*mean(fitdiag$g)),
                    zz+2*sqrt(fitdiag$s2_smooth*mean(fitdiag$tau2)*mean(fitdiag$g)),
                    zz-2*sqrt(fitdiag$s2*mean(fitdiag$tau2)),
                    zz+2*sqrt(fitdiag$s2*mean(fitdiag$tau2)),
                    zz-2*sqrt(1/precs_pred),
                    zz+2*sqrt(1/precs_pred), Yadj)), main = "fitdiag")
for (i in 1:nrow(Y))  lines(c(x), Y[i,]-y_avg)
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
# save(fitdiag, file = "rda/corr_err/fitdiag.rda")

xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(1/diag(covY)) ~ x) #precs for logl_g.R, 1/diag(covY) for logl_cov.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))


fitcov <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = nmcmc, Sigma_hat = covY/(16*sigma2_yavg))#, cov = "exp2")
# plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
plot(fitcov)
fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=6, precs_pred = precs_pred)

lmm <- lm(y_avg ~ x)
y_avg_e <- cbind(1,xx)%*%coefficients(lmm)

zz <- fitcov$mean-fitcov$mean

par(mfrow=c(1,1))
plot(fitcov$x_new, zz, type="l", col="blue", lwd=1.5, 
     ylim = range(c(zz-2*sqrt(fitcov$s2_smooth),
                    zz+2*sqrt(fitcov$s2_smooth),
                    zz-2*sqrt(fitcov$s2),
                    zz+2*sqrt(fitcov$s2),
                    zz-2*sqrt(1/precs_pred),
                    zz+2*sqrt(1/precs_pred), Yadj)), main = "fitcov")
for (i in 1:nrow(Y))  lines(c(x), Y[i,]-y_avg)
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
# save(fitcov, file = "rda/corr_err/fitcov.rda")

# UQ for truth S by Bayes Theorem:
# S | \bar{Y} \sim(m, C) where
# C^{-1} = (Sigma_s(x,x'))^{-1} + 16*Sigma_epsilon^{-}
# m = C^{-1} * 16*Sigma_epsilon^{-} \bar{Y}

# In particular, using the MCMC output for the unknown parameters in Sigma_s(x,x'), 
# we can use the distribution [S | \bar{Y}] above to simulate a sample from posterior of S. 

# fitcov$theta_y
# fitcov$tau2
# fitcov$g
