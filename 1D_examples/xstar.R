# Make a LHS sample that goes into xstar.dat
# In order to run CosmicEmu/*/emu.exe

library(lhs)
library(multiplex)

set.seed(1)

# Some quick values to try running emulator if necessary
# sample_params <- c(0.13, 0.022, 0.8, 0.7, 0.95, -1, -1, 0.005)

# Do you want conditional main effects (TRUE) or average main effects (FALSE)?
fix_0.5 <- T

G <- 11 # number of grid values to try: 0.0, 0.1, ..., 1
MC <- 50 # number of MC for each grid value

# Don't need Monte Carlo variation if conditional effects
if(fix_0.5) MC <- 1

n <- G*MC 

if(fix_0.5) xx <- matrix(0.5, n, 8) else  xx <- lhs::randomLHS(n, 8)

grid <- rep(seq(0, 1, length.out=G), each=MC)

xxx <- c()
for (i in 1:8) {
    xx_temp <- xx
    xx_temp[,i] <- grid
    xxx <- rbind(xxx, xx_temp)
}

# # If you want a 3^8 full factorial design (e.g.), 
# # use this code to make xxx instead!
# vals <- c(0, 0.5, 1)
# xxx <- expand.grid(vals, vals, vals, vals, vals, vals, vals, vals, 0)

omega_m <- xxx[,1]*(0.155-0.12)+0.12
omega_b <- xxx[,2]*(0.0235-0.0215)+0.0215
sigma_8 <- xxx[,3]*(0.9-0.7)+0.7
h <- xxx[,4]*(0.85-0.55)+0.55
n_s <- xxx[,5]*(1.05-0.85)+0.85
w_0 <- xxx[,6]*(-0.7--1.3)+-1.3
f_wa0 <- xxx[,7]*(1.2899999-0.3000001)+0.3000001 # f(w0, wa) = -(w_0+w_a)^(1/4)
w_a <- -f_wa0^4 - w_0
omega_v <- xxx[,8]*(0.01-0.0)+0.0
z <- 0

xstar <- cbind(omega_m, omega_b, sigma_8, h, n_s, w_0, w_a, omega_v, z)

# ensure all design points now make good cosmologies
w_a <- xstar[,"w_a"]
w_0 <- xstar[,"w_0"]
ind1 <- w_0 + w_a < -0.008099
ind2 <- (-w_0-w_a)^(1/4) < 1.29 & (-w_0-w_a)^(1/4) > 0.3
(inds <- which(!(ind1 & ind2)))

if(length(inds)>0) stop("some design points are still bad cosmologies")

write.dat(xstar, path = "~/CosmicEmu/2022-Mira-Titan-IV/P_tot/")

