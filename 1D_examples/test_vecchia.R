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
precs <- prec_lowres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
y <- log10(pk2[index_list$lowres.ix, bte])
y_avg <- rowMeans(log10(pk2[index_list$lowres.ix, 3:18]))
y_hi <- log10(pk2[index_list$lowres.ix, 19])

# save var(y_avg) to adjust prec values in the final fit call
sigma2_y <- var(y)
sigma2_yavg <- var(y_avg)

x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean(y))/sd(y)
y_avg <- (y_avg - mean(y_avg))/sd(y_avg)
y_hi <- (y_hi - mean(y_hi))/sd(y_hi)

#######################
# run for actual data #
#######################

# # original fit for the average of the 16 low res avg (no vectorized nugget)
# fit <- fit_two_layer(x, y_avg, cov = "matern", v=2.5, nmcmc = 50000, vecchia = T)
# save(fit,file=paste0("rda/fit3_avg_pm0_log10_vecc.rda"))
# fit <- fit_two_layer(x, y_avg, cov = "exp2", nmcmc = 50000, vecchia = T)
# save(fit,file=paste0("rda/fit4_avg_pm0_log10_vecc.rda"))

fit4 <- fit_two_layer_SW(x = x, y = y_avg, nmcmc = 252500, precs = (precs*sigma2_yavg*16), cov = "matern", v=2.5)
fit4 <- trim_SW(fit4, 2500, 2)
save(fit4, file = "rda/prec_vector/fitavg.rda")

fit4v <- fit_two_layer_SW(x, y_avg, nmcmc = 252500, precs = (precs*sigma2_yavg*16), vecchia = T)
fit4v <- trim_SW(fit4v, 2500, 2)
save(fit4v, file = "rda/prec_vector/fitavg_vec.rda")

fit4ve <- fit_two_layer(x, y_avg, nmcmc = 252500, vecchia = T)
fit4ve <- trim_SW(fit4ve, 2500, 2)
save(fit4ve, file = "rda/prec_vector/fitavg_vec_scalarg.rda")

fit4$time
fit4v$time
fit4ve$time
