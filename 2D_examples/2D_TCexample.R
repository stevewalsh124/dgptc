#################################
# Parameter estimation for TCs  #
# January 31, 2022              #
# Steve Walsh, Annie Sauer      #
#################################

library(raster)
library(deepgp) # version >= 0.3.0

tic <- proc.time()[3]

tc_files <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/", full.names = T)

# storm to evaluate (11 is smallest)
ste <- 11

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

tc <- read.csv(tc_files[ste], row.names = 1)

# scaled lon and lat to be in [0,1]
tc$xs <- (tc$x - min(tc$x))/(max(tc$x)-min(tc$x))
tc$ys <- (tc$y - min(tc$y))/(max(tc$y)-min(tc$y))

# set training and test sets
train <- sample(1:nrow(tc), 500)
test <- (1:nrow(tc))[-train]

# subset training data
tc_samp <- tc[train,]
x <- cbind(tc_samp$xs, tc_samp$ys)
y <- tc_samp$value
# plot(rasterFromXYZ(data.frame(cbind(x,y))))

# Locations for predictions
tc_pred <- tc[test,]
xx <- cbind(tc_pred$xs, tc_pred$ys)

niters <- 40000

# Fit two-layer DGP (exponential cov fn)
fit <- fit_two_layer(x, y, nmcmc = niters, cov = "matern", v=0.5, vecchia = T)
fit <- trim(fit, 1000, 1) # retain 2500 samples
fit <- predict(fit, xx)

# Fit two-layer DGP (Matern, v=5/2)
# fit2 <- fit_two_layer(x, y, nmcmc = niters, cov = "matern", v=2.5, vecchia = T)
# fit2 <- trim(fit2, 1000, 1)
# fit2 <- predict(fit2, xx)

# Combine results
pred <- data.frame(xx = xx, mean = fit$mean, s2 = fit$s2_smooth)
# pred2 <- data.frame(xx = xx, mean = fit2$mean, s2 = fit2$s2_smooth)
write.csv(pred, paste0("csv/",ste,".csv"), row.names = FALSE)

# get range for plots
rg <- range(c(tc$value, pred$mean)) #, pred2$mean
s2_rg <- range(c(pred$s2)) #, pred2$s2

# compare predictions via plots
pdf(paste0("expntl_vs_gaussian_cov_fn_vecchia_",niters,".pdf"))
par(mfrow=c(1,2))
plot(rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value)), zlim=rg)
plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$mean)), zlim=rg)
# plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$mean)), zlim=rg)

plot(rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value)))
plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$s2)), zlim=s2_rg)
# plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$s2)), zlim=s2_rg)

par(mfrow=c(1,3))
plot(fit$theta_w[,1], type="l")
plot(fit$theta_w[,2], type="l")
plot(fit$theta_y, type="l")

# plot(fit2$theta_w[,1], type="l")
# plot(fit2$theta_w[,2], type="l")
# plot(fit2$theta_y, type="l")
dev.off()

toc <- proc.time()[3]

toc - tic
