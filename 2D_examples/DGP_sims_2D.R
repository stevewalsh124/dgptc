##########################
# Simulating deep GPs
# Steve Walsh
# January 28 2022
##########################

# Code based on that of M.A.R. Ferreira, 2016.

pdf("pdf/deform_2D_sims.pdf")

library(fields) # image.plot
library(Morpho) # deformGrid2d
library(marmap) # griddify
library(raster) # plot(.raster)

# How many simulations do you want?
nsims <- 15
# # of pixels per simulation? (11^2, 21^2, 31^2, etc.)
npix <- 41^2

################################
# Gaussian covariance function #
################################

covfunc.Gaussian <- function(t,phi,sigma2) {sigma2 * exp(-t^2/phi^2)}

matrix.sqrt <- function(H)
{
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-10)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  

  H.sqrt
}

# With nugget effect
plot.gaus <- function(t, phi=1, sigma2=1, tau2=0){
  covmatrix <- diag(tau2,nrow(t)) + covfunc.Gaussian(t,phi,sigma2)
  x <- matrix.sqrt(covmatrix) %*% rnorm(nrow(s))
  z <- matrix(x,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
  image.plot(z, main = bquote(phi == .(phi) ~ ",  " ~ sigma^2 == 
                                .(sigma2) ~ ",  " ~ tau^2 == .(tau2)),
             cex.main=1.5)
}


####################
# Plot Simulations #
####################

xaxis <- seq(0,1, length.out = sqrt(npix))
yaxis <- seq(0,1, length.out = sqrt(npix))

# original lat lon coords
s <- matrix(NA,nrow=length(xaxis)*length(yaxis),ncol=2)
for (i in 1:length(xaxis)) for (j in 1:length(yaxis)) s[i+(j-1)*length(xaxis),] <- c(xaxis[i],yaxis[j])
# plot(s)

# first layer...
t <- as.matrix(dist(s))
ss <- rowSums(s)

for (i in 1:nsims) {
  
  print(i)
  
  par(mfrow=c(2,2))
  # Currently the same for each layer; can change for individual layers
  phi_y <- 1/3
  phi_w <- 0.01*i^3
  sigma2 <- 1
  tau2 <- .Machine$double.eps^0.5
  
  # ... first layer continued
  covmatrix <- diag(tau2,nrow(s)) + covfunc.Gaussian(t,phi_w,sigma2)
  w1 <- t(chol(covmatrix)) %*% rnorm(nrow(s)) + s[,1] # + ss
  w2 <- t(chol(covmatrix)) %*% rnorm(nrow(s)) + s[,2] # + ss
  w <- cbind(w1,w2)
  
  # second layer
  t2 <- as.matrix(dist(w))
  covmatrix2 <- diag(tau2,nrow(s)) + covfunc.Gaussian(t2,phi_y,sigma2)
  y <- t(chol(covmatrix2)) %*% rnorm(nrow(s))
  
  # # third layer
  # t3 <- as.matrix(dist(w))
  # covmatrix3 <- diag(tau2,nrow(s)) + covfunc.Gaussian(t3,phi,sigma2)
  # z <- t(chol(covmatrix2)) %*% rnorm(nrow(w))
  
  # output; change "matrix(w" to "matrix(z" if using third layer
  yp <- matrix(y,nrow=length(xaxis),ncol=length(yaxis),byrow=FALSE)
  nug <- rnorm(prod(dim(y)), sd = sqrt(tau2))
  
  # establish which points in x (and y) are reference points
  # which low (wl) is (0.5,0.1) and which high (wh) is (0.5,0.9)
  wl <- which(abs(s[,1]-0.5) < 1e-3 & abs(s[,2]-0.1) < 1e-3)
  wh <- which(abs(s[,1]-0.5) < 1e-3 & abs(s[,2]-0.9) < 1e-3)
  
  # If reference point pair are not in training set, skip
  if(length(wl)==0 | length(wh)==0) next

  zz <- s[wl,]
  zo <- s[wh,]
  # points(zz[1],zz[2], lwd=3) # add to plot above
  # points(zo[1],zo[2], lwd=3) # add to plot above
  
  # this finds the length of the original reference vector
  zzo <- zz - zo
  ref.scale <- sqrt(t(zzo)%*%zzo)

  # Don't plot: the original (unadjusted) deformation
  # deformGrid2d(s, w, ngrid=25, pch=19, main=paste(i, "old"), gridcol = "black")
  
  # Take these deformed points (y) and rotate/scale/translate them
  # to match orig lat/lon points (x) as well as possible
  st <- t(s)
  wt <- t(w)
  
  # First: undo rotation
  ang <- function(u,v){acos((t(u)%*%v)/((sqrt(t(u)%*%u) * sqrt(t(v)%*%v))))}
  u <- st[,wh] - st[,wl]
  v <- wt[,wh] - wt[,wl]
  phi <- ang(u,v)
  rotMtx <- matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
  y.rot <- rotMtx %*% wt
  
  # Second: undo scale
  # find the length of the deformed reference vector
  y.ref <- wt[,wl] - wt[,wh]
  scale.y <- sqrt(t(y.ref)%*%y.ref)
  # numerator comes from length of vector with points zz, zo
  scale.back =  c(ref.scale/scale.y)#1/sqrt(sum(y.translate[,2]^2))
  y.rot.scale = y.rot * scale.back
  
  # Third: Translate back
  translation <- y.rot.scale[,wl] - st[,wl]
  y.rot.scale.tran <- y.rot.scale - translation
  ww <- t(y.rot.scale.tran)

  # show plots of deformation
  # same as plot(rasterFromXYZ(cbind(s,w1)))
  image.plot(matrix(w1,length(xaxis),length(yaxis)), main = "x to w1") 
  image.plot(matrix(w2,length(xaxis),length(yaxis)), main = "x to w2")
  
  # plot the complete deformation of geographic points
  deformGrid2d(s, ww, ngrid=25, pch=19, main=paste("X to W"), gridcol = "black", lines = F)
  
  # Deep GP is made by iso GP of (iso GP of locns) = iso GP of (deformed space)
  irreg <- as.data.frame(cbind(ww,y))
  colnames(irreg) <- c("lon","lat","y")
  reg <- griddify(irreg, nlon = 40, nlat = 60)
  plot(reg, main = "W to Y")
  
  par(mfrow=c(1,2))
  # compare the stationary, isotropic GP on the deformed space...
  plot(reg, main = "W to Y (iso, stat)")
  # ... with the corresponding nonstationary GP on geographic points
  plot(rasterFromXYZ(cbind(s,y)), main = "X to Y (nonst)")
  # # same as above
  # image.plot(yp + matrix(nug, length(xaxis),length(yaxis)),
  #            main = bquote("X to Y\n" ~ phi == .(round(phi,1)) ~ ",  " ~ sigma^2 == .(round(sigma2,1)) ~ ",  "
  #                          ~ tau^2 == .(round(tau2,4))), cex.main=1.5)
}

dev.off()
