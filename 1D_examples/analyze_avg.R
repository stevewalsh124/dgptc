###############################
# analyze average low res run #
# 1D cosmology example        #
# April 28 2022                #
###############################

library(deepgp)

nburn <- 0
kth <- 1 # every kth sample

wl <- 1; wh <- 161
ref.scale <- 1

pdf(paste0("pdf/average_runs_",kth,".pdf"))

# can run this chunk with fit, fit2, fit3, and fit4
load("rda/prec_vector/fitavg_thin10.rda")
fit <- fit4

plot(fit)

fit <- trim(fit, nburn, kth)
x <- fit$x

for (j in 1:fit$nmcmc) {
  y <- c(fit$w[[j]])
  
  # Second: undo scale
  # find the length of the deformed reference vector
  y.ref <- y[wl] - y[wh]
  scale.y <- sqrt(t(y.ref)%*%y.ref)
  # numerator comes from length of vector with points zz, zo
  scale.back =  c(ref.scale/scale.y)#1/sqrt(sum(y.translate[,2]^2))
  y.rot.scale = y * scale.back
  
  # Third: Translate back
  translation <- y.rot.scale[wl] - x[wl]
  y.rot.scale.tran <- y.rot.scale - translation
  wstar <- y.rot.scale.tran
  
  # flip the graph if necessary in either direction
  # remove the reflections
  if(wstar[wl]>wstar[wh]) {wstar <- (wstar-wstar[wl])/(wstar[wh]-wstar[wl])}
  fit$w[[j]] <- wstar
}

ws <- do.call(rbind, fit$w)
w_avg <- colMeans(ws)
w_lb <- apply(ws, 2, function(x) quantile(x, 0.025))
w_ub <- apply(ws, 2, function(x) quantile(x, 0.975))

par(mfrow=c(1,1))
plot(x,fit$w[[j]], type="l", main = "warpings, with mean and 95% CI")
for(j in 1:fit$nmcmc) lines(x, fit$w[[j]])
lines(x, w_avg, col="blue")
lines(x, w_lb, col="blue")
lines(x, w_ub, col="blue")
abline(0, 1, col="red")

par(mfrow=c(2,2))
plot(fit$theta_w, type="l", main = paste("theta_w", mean(fit$theta_w)))
plot(fit$theta_y, type="l", main = paste("theta_y", mean(fit$theta_y)))
plot(fit$tau2, type="l", main = paste("tau2", mean(fit$tau2)))
plot(fit$g, type="l", main = paste("g", mean(fit$g)))

hist(fit$theta_w, main = paste("theta_w", mean(fit$theta_w)))
hist(fit$theta_y, main = paste("theta_y", mean(fit$theta_y)))
hist(fit$tau2, main = paste("tau2", mean(fit$tau2)))
hist(fit$g, main = paste("g", mean(fit$g)))

plot(density(fit$theta_w), main = paste("theta_w", mean(fit$theta_w)))
plot(density(fit$theta_y), main = paste("theta_y", mean(fit$theta_y)))
plot(density(fit$tau2), main = paste("tau2", mean(fit$tau2)))
plot(density(fit$g), main = paste("g", mean(fit$g)))

dev.off()
