# Make a LHS sample that goes into xstar.dat
# In order to run CosmicEmu/*/emu.exe

library(lhs)
library(multiplex)

set.seed(1)

# Some quick values to try running emulator if necessary
# sample_params <- c(0.13, 0.022, 0.8, 0.7, 0.95, -1, -1, 0.005)

G <- 11 # number of grid values to try: 0.0, 0.1, ..., 1
MC <- 20 # number of MC for each grid value

#!#!# this is hacky; trying to get joint G*MC 
# after removing 
n <- G*(MC*1.15) 
xx <- lhs::randomLHS(n, 8)

# rm design values outside of joint range for w_0 and w_a
w_0 <- xx[,6]*(-0.7--1.3)+-1.3
w_a <- xx[,7]*(1.28--1.73)+-1.73
ind1 <- w_0 + w_a < -0.0081
ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
xx <- xx[ind1 & ind2,]  #!#!#
if(nrow(xx) < G*MC) stop("do more runs...")
xx <- xx[1:(G*MC),]

grid <- rep(seq(0, 1, length.out=G), each=MC)

xxx <- c()
for (i in 1:8) {
    xx_temp <- xx
    xx_temp[,i] <- grid
    xxx <- rbind(xxx, xx_temp)
}

omega_m <- xxx[,1]*(0.155-0.12)+0.12
omega_b <- xxx[,2]*(0.0235-0.0215)+0.0215
sigma_8 <- xxx[,3]*(0.9-0.7)+0.7
h <- xxx[,4]*(0.85-0.55)+0.55
n_s <- xxx[,5]*(1.05-0.85)+0.85
w_0 <- xxx[,6]*(-0.7--1.3)+-1.3
w_a <- xxx[,7]*(1.28--1.73)+-1.73
omega_v <- xxx[,8]*(0.01-0.0)+0.0
z <- 0

xstar <- cbind(omega_m, omega_b, sigma_8, h, n_s, w_0, w_a, omega_v, z)

# find which design points will make bad cosmologies
ind1 <- w_0 + w_a < -0.0081
ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
inds <- which(!(ind1 & ind2))

# randomly replace w_0 or w_a but keep griddedness
for (ff in 1:length(inds)) {
  fff <- inds[ff]
  w_0 <- xstar[fff,"w_0"]
  w_a <- xstar[fff,"w_a"]
  ind1 <- w_0 + w_a < -0.0081
  ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
  if(fff > 1100 & fff <= 1320){
    # fix w_0, change w_a
    while(!ind1 | !ind2){
      w_a <- runif(1)*(1.28--1.73)+-1.73
      ind1 <- w_0 + w_a < -0.0081
      ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
    }
    xstar[fff,"w_a"] <- w_a
    
  } else {
    
    if(fff > 1320 & fff <= 1540){
      while (!ind1 | !ind2) {
        # fix w_a, change w_0
        w_0 <- runif(1)*(-0.7--1.3)+-1.3
        ind1 <- w_0 + w_a < -0.0081
        ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
      }
      xstar[fff,"w_0"] <- w_0
      
    } else {
      stop("something is wrong, should only be changing w_0 or w_a")
    }
  }
}

# ensure all design points now make good cosmologies
w_a <- xstar[,"w_a"]
w_0 <- xstar[,"w_0"]
ind1 <- w_0 + w_a < -0.0081
ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
inds <- which(!(ind1 & ind2))

if(length(inds)>0) stop("some design points are still bad cosmologies")

write.dat(xstar, path = "~/CosmicEmu/2022-Mira-Titan-IV/P_tot/")

