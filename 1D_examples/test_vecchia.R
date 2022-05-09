source("logl_cov.R")

nmcmc <- 25000
nburn <- 5000
kth <- 2

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

step <- 499
i <- 1 # Model 1, choose from 000-111
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

# obstain sample cov mat after rescaling
covY <- cov(Y)

# get predictions for xx and have corresponding prec info
xx <- setdiff(seq(0,1,by=.01), x)
lmfit <- lm(log10(1/diag(covY)) ~ x) #precs for logl_g.R
betahat <- coef(lmfit)
precs_pred <- as.numeric(10^(cbind(1,xx) %*% betahat))

#######################
# run for actual data #
#######################

# fit4 <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = nmcmc, precs = precs*sigma2_yavg*16)
# fit4 <- trim_SW(fit4, nburn, kth)
# fit4 <- predict.dgp2_SW(fit4, xx, precs_pred = (precs_pred*sigma2_yavg*16), cores=2)
# plot(fit4)
# plot(fit4$x_new, fit4$mean, type="l", col="blue")
# lines(fit4$x_new, fit4$mean-30*sqrt(fit4$s2_smooth))
# lines(fit4$x_new, fit4$mean+30*sqrt(fit4$s2_smooth))
# lines(fit4$x_new, fit4$mean-30*sqrt(fit4$s2))
# lines(fit4$x_new, fit4$mean+30*sqrt(fit4$s2))
# save(fit4, file = "rda/prec_vector/new/fit4mat.rda")

fitcov <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = nmcmc, Sigma_hat = covY)#, cov = "exp2")
plot(fitcov)
fitcov <- trim_SW(fitcov, nburn, kth)
plot(fitcov)
fitcov <- predict.dgp2_SW(object = fitcov, xx, cores=6, precs_pred = precs_pred)

par(mfrow=c(1,1))
plot(fitcov$x_new, fitcov$mean, type="l", col="blue", lwd=1.5)
for (i in 1:nrow(Y))  lines(c(x), Y[i,])
lines(fitcov$x_new, fitcov$mean, col="blue", lwd=1.5)
lines(fitcov$x_new, fitcov$mean-2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, fitcov$mean+2*sqrt(fitcov$s2_smooth), col=2, lwd=1.5)
lines(fitcov$x_new, fitcov$mean-2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, fitcov$mean+2*sqrt(fitcov$s2), col=2, lwd=1.5)
lines(fitcov$x_new, fitcov$mean-2*sqrt(1/precs_pred), col=3, lwd=1.5)
lines(fitcov$x_new, fitcov$mean+2*sqrt(1/precs_pred), col=3, lwd=1.5)
legend("topright", col=c("blue",2:3), legend = c("mean","UQ", "precs (low res)"), lty=1, lwd=1.5)
# save(fitcov, file = "rda/prec_vector/new/fitcovmatvecc.rda")
