source("logl_g.R")

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
vars <- 1/prec_lowres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
y <- log10(pk2[index_list$lowres.ix, bte])
y_avg <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]))
y_hi <- log10(pk2[index_list$lowres.ix, 19])

# save var(y_avg) to adjust prec values in the final fit call
sigma2_y <- var(y_avg)

x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean(y))/sd(y)
y_avg <- (y_avg - mean(y_avg))/sd(y_avg)
y_hi <- (y_hi - mean(y_hi))/sd(y_hi)

#############################
# pdf: eg of different runs #
#############################

pdf("pdf/low_hi_avg_runs_log10.pdf")
plot(x,y, type="l")
for (j in 4:18) {
  x <- log10(k[index_list$lowres.ix])
  y <- log10(pk2[index_list$lowres.ix, j])
  
  # standardize inputs and outputs
  x <- (x-min(x))/(max(x)-min(x))
  y <- (y - mean(y))/sd(y)
  lines(x,y)
}
lines(x,y_avg, col="red", lwd=2)
lines(x,y_hi, col="lightblue", lwd=2)
legend("topright",c("hi","avg"), col=c("lightblue","red"), lty=1)
dev.off()

################################
# Check equality for full lkhd #
################################

# ensure results are the same for scalar g
set.seed(1)
fit <- fit_two_layer(x, y, nmcmc = 10, true_g = 1e-4)
set.seed(1)
fit2 <- fit_two_layer_SW(x, y, nmcmc = 10, true_g = 1e-4) #doesn't actually use SW version (scalar)
set.seed(1)
fit3 <- fit_two_layer_SW(x, y, nmcmc = 10, true_g = rep(1e-4, nrow(as.matrix(x))))

all.equal(fit$theta_y, fit2$theta_y)
all.equal(fit$theta_w, fit2$theta_w)
all.equal(fit$tau2, fit2$tau2)
all.equal(fit$g, fit2$g)

all.equal(fit$theta_y, fit3$theta_y)
all.equal(fit$theta_w, fit3$theta_w)
all.equal(fit$tau2, fit3$tau2)
for(i in 1:nrow(as.matrix(x))) if(!all.equal(fit$g, fit3$g[,i])) stop("g's not equal :(")

########################
# Check Vecchia approx #
########################

set.seed(1)
fit <- fit_two_layer(x, y, nmcmc = 10, true_g = 1e-4, vecchia = T)
set.seed(1)
fit2 <- fit_two_layer_SW(x, y, nmcmc = 10, true_g = 1e-4, vecchia = T) #doesn't actually use SW version (scalar)
set.seed(1)
fit3 <- fit_two_layer_SW(x, y, nmcmc = 10, true_g = rep(1e-4, nrow(as.matrix(x))), vecchia = T)

all.equal(fit$theta_y, fit2$theta_y)
all.equal(fit$theta_w, fit2$theta_w)
all.equal(fit$tau2, fit2$tau2)
all.equal(fit$g, fit2$g)

all.equal(fit$theta_y, fit3$theta_y)
all.equal(fit$theta_w, fit3$theta_w)
all.equal(fit$tau2, fit3$tau2)
for(i in 1:nrow(as.matrix(x))) if(!all.equal(fit$g, fit3$g[,i])) stop("g's not equal :(")

#######################
# run for actual data #
#######################

# # original fit for the average of the 16 low res avg (no vectorized nugget)
# fit <- fit_two_layer(x, y_avg, cov = "matern", v=2.5, nmcmc = 50000, vecchia = T)
# save(fit,file=paste0("rda/fit3_avg_pm0_log10_vecc.rda"))
# fit <- fit_two_layer(x, y_avg, cov = "exp2", nmcmc = 50000, vecchia = T)
# save(fit,file=paste0("rda/fit4_avg_pm0_log10_vecc.rda"))

fit4 <- fit_two_layer_SW(x, y_avg, nmcmc = 51000, true_g = vars*sigma2_y/16)
fit4 <- trim_SW(fit4, 1000)
save(fit4, file = "rda/g_vector/fit4.rda")

fit5 <- fit_two_layer_SW(x, 3*y_avg, nmcmc = 51000, true_g = vars*sigma2_y/16)
fit5 <- trim_SW(fit5, 1000)
save(fit5, file = "rda/g_vector/fit5.rda")

fit6 <- fit_two_layer_SW(x, y_avg/3, nmcmc = 51000, true_g = vars*sigma2_y/16)
fit6 <- trim_SW(fit6, 1000)
save(fit6, file = "rda/g_vector/fit6.rda")

fit4v <- fit_two_layer_SW(x, y_avg, nmcmc = 51000, true_g = vars*sigma2_y/16, vecchia = T)
fit4v <- trim_SW(fit4v, 1000)
save(fit4v, file = "rda/g_vector/fit4v.rda")

fit5v <- fit_two_layer_SW(x, 3*y_avg, nmcmc = 51000, true_g = vars*sigma2_y/16, vecchia = T)
fit5v <- trim_SW(fit5v, 1000)
save(fit5v, file = "rda/g_vector/fit5v.rda")

fit6v <- fit_two_layer_SW(x, y_avg/3, nmcmc = 51000, true_g = vars*sigma2_y/16, vecchia = T)
fit6v <- trim_SW(fit6v, 1000)
save(fit6v, file = "rda/g_vector/fit6v.rda")

fit4ve <- fit_two_layer(x, y_avg, nmcmc = 51000, true_g = vars*sigma2_y/16, vecchia = T)
fit4ve <- trim_SW(fit4ve, 1000)
save(fit4ve, file = "rda/g_vector/fit4vecchia_nonSW.rda")

fit5ve <- fit_two_layer(x, y_avg*3, nmcmc = 51000, true_g = vars*sigma2_y/16, vecchia = T)
fit5ve <- trim_SW(fit5ve, 1000)
save(fit5ve, file = "rda/g_vector/fit5vecchia_nonSW.rda")

fit6vvFAST <- fit_two_layer(x, y_avg/3, nmcmc = 5100, true_g = vars*sigma2_y/16, vecchia = T, ncores=16)
fit6ve <- trim_SW(fit6ve, 1000)
save(fit6ve, file = "rda/g_vector/fit6vecchia_nonSW.rda")



c(mean(fit4$theta_w), mean(fit5$theta_w), mean(fit6$theta_w),
  mean(fit4v$theta_w), mean(fit5v$theta_w), mean(fit6v$theta_w),
  mean(fit4ve$theta_w), mean(fit5ve$theta_w), mean(fit6ve$theta_w))

boxplot(c(fit4$theta_w), c(fit5$theta_w), c(fit6$theta_w),
        c(fit4v$theta_w), c(fit5v$theta_w), c(fit6v$theta_w),
        c(fit4ve$theta_w), c(fit5ve$theta_w), c(fit6ve$theta_w), main = "theta_w")

c(mean(fit4$tau2), mean(fit5$tau2), mean(fit6$tau2),
  mean(fit4v$tau2), mean(fit5v$tau2), mean(fit6v$tau2),
  mean(fit4ve$tau2), mean(fit5ve$tau2), mean(fit6ve$tau2))

boxplot(fit4$tau2, fit5$tau2, fit6$tau2,
        fit4v$tau2, fit5v$tau2, fit6v$tau2,
        fit4ve$tau2, fit5ve$tau2, fit6ve$tau2, main = "tau2", log="y")

c(mean(fit4$theta_y), mean(fit5$theta_y), mean(fit6$theta_y),
  mean(fit4v$theta_y), mean(fit5v$theta_y), mean(fit6v$theta_y),
  mean(fit4ve$theta_y), mean(fit5ve$theta_y), mean(fit6ve$theta_y))

boxplot(fit4$theta_y, fit5$theta_y, fit6$theta_y,
        fit4v$theta_y, fit5v$theta_y, fit6v$theta_y,
        fit4ve$theta_y, fit5ve$theta_y, fit6ve$theta_y, main = "theta_y", log="y")

c(mean(fit4$tau2/fit4$theta_y), mean(fit5$tau2/fit5$theta_y), mean(fit6$tau2/fit6$theta_y),
  mean(fit4v$tau2/fit4v$theta_y), mean(fit5v$tau2/fit5v$theta_y), mean(fit6v$tau2/fit6v$theta_y),
  mean(fit4ve$tau2/fit4ve$theta_y), mean(fit5ve$tau2/fit5ve$theta_y), mean(fit6ve$tau2/fit6ve$theta_y))

boxplot(fit4$tau2/fit4$theta_y, fit5$tau2/fit5$theta_y, fit6$tau2/fit6$theta_y,
        fit4v$tau2/fit4v$theta_y, fit5v$tau2/fit5v$theta_y, fit6v$tau2/fit6v$theta_y,
        fit4ve$tau2/fit4ve$theta_y, fit5ve$tau2/fit5ve$theta_y, fit6ve$tau2/fit6ve$theta_y, 
        main = "tau2/theta_y", log="y")

c(mean(fit4$time), mean(fit5$time), mean(fit6$time),
  mean(fit4v$time), mean(fit5v$time), mean(fit6v$time),
  mean(fit4ve$time), mean(fit5ve$time), mean(fit6ve$time))

boxplot(fit4$time, fit5$time, fit6$time,
        fit4v$time, fit5v$time, fit6v$time,
        fit4ve$time, fit5ve$time, fit6ve$time, main = "time")

c(mean(fit4$g), mean(fit5$g), mean(fit6$g),
  mean(fit4v$g), mean(fit5v$g), mean(fit6v$g),
  mean(fit4ve$g), mean(fit5ve$g), mean(fit6ve$g))

boxplot(fit4$g[1,], fit5$g[1,], fit6$g[1,],
        fit4v$g[1,], fit5v$g[1,], fit6v$g[1,],
        fit4ve$g, fit5ve$g, fit6ve$g, main = "g")