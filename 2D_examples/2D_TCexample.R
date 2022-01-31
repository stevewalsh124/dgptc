
library(raster)
library(deepgp) # version >= 0.3.0

tic <- proc.time()[3]

tc <- read.csv("../data/errordf_2005katrina_700deg.csv", row.names = 1)

# create training data
tc_samp <- tc[sample(1:nrow(tc), 500),]
x <- cbind(tc_samp$x, tc_samp$y)
y <- tc_samp$value
# plot(rasterFromXYZ(data.frame(cbind(x,y))))

# Locations for predictions
xx <- cbind(tc$x, tc$y)

niters <- 10000

# Fit two-layer DGP (exponential cov fn)
fit <- fit_two_layer(x, y, nmcmc = niters, cov = "matern", v=0.5)
fit <- trim(fit, 1000, 1) # retain 2500 samples
fit <- predict(fit, xx)

# Fit two-layer DGP (Gaussian cov fn)
fit2 <- fit_two_layer(x, y, nmcmc = niters, cov = "exp2")
fit2 <- trim(fit2, 1000, 1) # retain 2500 samples
fit2 <- predict(fit2, xx)

# Combine results
pred <- data.frame(xx = xx, mean = fit$mean, s2 = fit$s2_smooth)
pred2 <- data.frame(xx = xx, mean = fit2$mean, s2 = fit2$s2_smooth)
# write.csv(pred, "piecewise_predictions.csv", row.names = FALSE)

# get range for plots
rg <- range(c(tc$value, pred$mean, pred2$mean))
s2_rg <- range(c(pred$s2, pred2$s2))

# compare predictions via plots
pdf(paste0("expntl_vs_gaussian_cov_fn_",niters,".pdf"))
par(mfrow=c(1,3))
plot(rasterFromXYZ(tc[,c(2,3,1)]), zlim=rg)
plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$mean)), zlim=rg)
plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$mean)), zlim=rg)

plot(rasterFromXYZ(tc[,c(2,3,1)]))
plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$s2)), zlim=s2_rg)
plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$s2)), zlim=s2_rg)
dev.off()

toc <- proc.time()[3]

toc - tic
