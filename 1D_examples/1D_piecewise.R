
library(deepgp) # version >= 0.3.0

piecewise <- function(x) {
  if (x > -0.5 & x < 0) return(exp(4 + 5 / (2 * x * (10 * x + 5))))
  if (x >= 0.2 & x <= 0.3) return(1)
  if (x > 0.3 & x <= 0.4) return(-1)
  return(0)
}

n <- 300
sig2 <- 0.01
x <- runif(n, -0.5, 0.5)
y <- sapply(x, piecewise) + rnorm(n, 0, sqrt(sig2))

m <- 1001
xx <- seq(-0.5, 0.5, length = m)
yy <- sapply(xx, piecewise)
# plot(x, y); lines(xx, yy, col = 2)

# Fit two-layer DGP
fit <- fit_two_layer(x, y, nmcmc = 20000, cov = "exp2")
fit <- trim(fit, 15000, 2) # retain 2500 samples
fit <- predict(fit, xx)

# Combine results
pred <- data.frame(xx = xx, yy = yy, mean = fit$mean, s2 = fit$s2_smooth)
write.csv(pred, "piecewise_predictions.csv", row.names = FALSE)

