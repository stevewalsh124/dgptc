
source("metrics.R") # requires "verification" package
shade_color <- rgb(col2rgb("blue")[1], col2rgb("blue")[2], col2rgb("blue")[3], 
                   max = 255, alpha = 45)

step <- read.csv("step_predictions.csv")
piecewise <- read.csv("piecewise_predictions.csv")

# Evaluate prediction metrics
results <- data.frame(
  step = c(MAPE(step$yy, step$mean),
           RMSPE(step$yy, step$mean),
           CRPS(step$yy, step$mean, step$s2),
           IS95(step$yy, step$mean, step$s2)),
  piecewise = c(MAPE(piecewise$yy, piecewise$mean),
                RMSPE(piecewise$yy, piecewise$mean),
                CRPS(piecewise$yy, piecewise$mean, piecewise$s2),
                IS95(piecewise$yy, piecewise$mean, piecewise$s2)))
rownames(results) <- c("MAPE", "RMSPE", "CRPS", "IS95")
print(results)

# Plot step model
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
lower <- step$mean - 2 * sqrt(step$s2)
upper <- step$mean + 2 * sqrt(step$s2)
plot(step$xx, step$yy, type = "l", xlab = "x", ylab = "y",
     ylim = c(-1.5, 1.5), main = "Y1")
polygon(c(step$xx, rev(step$xx)), c(lower, rev(upper)), 
        col = shade_color, border = NA)
lines(step$xx, step$mean, col = "blue")

# Plot piecewise model
lower <- piecewise$mean - 2 * sqrt(piecewise$s2)
upper <- piecewise$mean + 2 * sqrt(piecewise$s2)
plot(piecewise$xx, piecewise$yy, type = "l", xlab = "x", ylab = "y",
     ylim = c(-1.5, 1.5), main = "Y2")
polygon(c(piecewise$xx, rev(piecewise$xx)), c(lower, rev(upper)), 
        col = shade_color, border = NA)
lines(piecewise$xx, piecewise$mean, col = "blue")
