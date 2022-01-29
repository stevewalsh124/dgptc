
# Copied from supplementary material for "Deep Compositional Spatial Models"
library(verification)

# Mean absolute prediction error
MAPE <- function(true, pred) 
  mean(abs(true - pred))

# Root mean squared prediction error
RMSPE <- function(true, pred) 
  sqrt(mean((true - pred)^2))

# Continuous ranked probability score, requires point-wise variances
CRPS <- function(true, pred, pred_var) 
  verification::crps(true, cbind(pred, sqrt(pred_var)))$CRPS

# Interval score of 95% prediction interval
# Revised to accept point-wise variance as input
IS95 <- function(true, pred, pred_var) {
  alpha <- 0.05
  pred95l <- pred - 2 * sqrt(pred_var)
  pred95u <- pred + 2 * sqrt(pred_var)
  ISs <- (pred95u - pred95l) + 2/alpha * (pred95l - true) * (true < pred95l) +
    2/alpha * (true - pred95u) * (true > pred95u)
  return(mean(ISs))
}
