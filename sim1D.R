###################################################
# Generate 1D data for deep GPs
# Author: Steve Walsh
# Date: January 27 2022
###################################################

x <- seq(-2,2,by=0.01)

# Different simulated functions
f1 <- function(t){cos(2*pi*3*t)*exp(-pi*t^2)}
f2 <- function(t){cos(2*pi*2*t-2*t^4)*exp(-pi*(t^2+t-0.5)^2)}
f3 <- function(t){sin(2*pi*3*t-2*t^4)*exp(-pi*(t-1)^2)}
f4 <- function(t){cos(2*pi*2*t-2*t^4)*exp(-pi*(t-0.5)^2)}

# Generate output
y1 <- f1(x)
y2 <- f2(x)
y3 <- f3(x)
y4 <- f4(x)

# Plot each function
par(mfrow=c(2,2))
plot(x,y1,type="l")
plot(x,y2,type="l")
plot(x,y3,type="l")
plot(x,y4,type="l")

# Combine data and export
data_1D <- cbind(x,y1,y2,y3,y4)
save(data_1D, file = "data/sims_1D.rda")
