#############################################
# Make rotation-invariant deformation plots #
# Steve Walsh                               #
# February 14 2022                          #
#############################################

tic <- proc.time()[3]

library(Morpho) #deformGrid2d
library(marmap) #griddify
library(raster) # rasterFromXYZ

niters <- 50058
k <- 10000 # every kth warping to plot
tp <- floor(niters/k) # total plots to produce (excluding original)

FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]

pdf(paste0("pdf/deform/deform_niters",niters,"_k",k,"_tp",tp,".pdf"))

load("rda/all_FL_tracks.rda")
load("rda/FL_rgs.rda")
for (i in fl) {
  all_FL_tracks[[i]][,1] <- (all_FL_tracks[[i]][,1] - FL_rgs[1]) / 
    (FL_rgs[2] - FL_rgs[1])
  all_FL_tracks[[i]][,2] <- (all_FL_tracks[[i]][,2] - FL_rgs[3]) / 
    (FL_rgs[4] - FL_rgs[3])
}

# Get reference points and scale
load("rda/FL_ref.rda")
zz <- FL_ref[1,]
zo <- FL_ref[2,]
sw <- FL_ref[3,]
se <- FL_ref[4,]
# points(zz[1],zz[2], lwd=3) # add to plot above
# points(zo[1],zo[2], lwd=3) # add to plot above

# this finds the length of the original reference vector
zzo <- zz - zo
ref.scale <- sqrt(t(zzo)%*%zzo)

# Loop over TCs
for (ste in 1:18) {
  load(paste0("rda/FL_fits/storm",ste,"_niters",niters,".rda")) #loads fit
  
  x <- t(fit$x)

  # The original non-invariant plots; shown in left side for comparison
  # par(mfrow=c(1,1))
  # deformGrid2d(fit$x, W[[1]], ngrid=25, pch=19, main=paste(ste, 1), gridcol = "black")
  
  # establish which points in x (and y) are reference points
  # which low (wl) is (0,0) and which high (wh) is (0,1)
  wl <- which(abs(x[1,]-FL_ref[1,1]) < 1e-3 & abs(x[2,]-FL_ref[1,2]) < 1e-3) # (0,0) point
  wh <- which(abs(x[1,]-FL_ref[2,1]) < 1e-3 & abs(x[2,]-FL_ref[2,2]) < 1e-3) # (0,1) point
  ww <- which(abs(x[1,]-sw[1]) < 1e-3 & abs(x[2,]-sw[2]) < 1e-3) # (0,1) point
  we <- which(abs(x[1,]-se[1]) < 1e-3 & abs(x[2,]-se[2]) < 1e-3) # (0,1) point
  
  # If reference point pair are not in training set, skip
  if(length(wl)==0 | length(wh)==0 | length(ww)==0 | length(we)==0) next
  print(ste)
  
  # For each simulated hidden layer: 
  for(i in k*(1:tp)){
    
    par(mfrow=c(2,2))
    # On the left: the original (unchanged) deformation
    deformGrid2d(fit$x, fit$w[[i]], ngrid=25, pch=19, main=paste(ste, i, "old"), gridcol = "black")
    
    # Take these deformed points (y) and rotate/scale/translate them
    # to match orig lat/lon points (x) as well as possible
    y <- t(fit$w[[i]])

    # First: undo rotation
    ang <- function(u,v){acos((t(u)%*%v)/((sqrt(t(u)%*%u) * sqrt(t(v)%*%v))))}
    u <- x[,wh] - x[,wl]
    v <- y[,wh] - y[,wl]
    phi <- ang(u,v)
    rotMtx <- matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
    y.rot <- rotMtx %*% y
    
    # Second: undo scale
    # find the length of the deformed reference vector
    y.ref <- y[,wl] - y[,wh]
    scale.y <- sqrt(t(y.ref)%*%y.ref)
    # numerator comes from length of vector with points zz, zo
    scale.back =  c(ref.scale/scale.y)#1/sqrt(sum(y.translate[,2]^2))
    y.rot.scale = y.rot * scale.back

    # Third: Translate back
    translation <- y.rot.scale[,wl] - x[,wl]
    y.rot.scale.tran <- y.rot.scale - translation
    wstar <- t(y.rot.scale.tran)
    
    # flip the graph if necessary in either direction
    if(wstar[ww,1]>wstar[we,1]){wstar[,1] <- (wstar[,1]-wstar[ww,1])/(wstar[we,1]-wstar[ww,1])}
    if(wstar[wh,2]>wstar[wl,2]){wstar[,2] <- (wstar[,2]-wstar[wl,2])/(wstar[wh,2]-wstar[wl,2])}
    
    # On right: plot a deformation that is rotation-invariant wrt original lat/lon
    deformGrid2d(fit$x, wstar, ngrid=25, pch=19, main=paste(ste, i), gridcol = "black")
    
    # Move the deformed points onto a regular grid
    yy <- fit$y
    irreg <- as.data.frame(cbind(wstar,yy))
    colnames(irreg) <- c("lon","lat","y")
    reg <- griddify(irreg, nlon = 40, nlat = 60)
    
    # Plot the deformed error field along with original points
    plot(fit$x[,1], fit$x[,2], pch = ".", cex = 0.3, asp=1,
         xlim = c(range(fit$x[,1], irreg$lon)), ylim = c(range(fit$x[,2], irreg$lat)))
    plot(reg, add=T, alpha=1)
    points(irreg$lon, irreg$lat, pch = ".", cex = 0.3, col = col2alpha(1, alpha = 0.1))
    
    # For reference, plot the original error field with track 
    plot(rasterFromXYZ( cbind(fit$x, fit$y)))
    points(all_FL_tracks[[fl[ste]]][,1], all_FL_tracks[[fl[ste]]][,2], pch=".")
    
  }
}

dev.off()

toc <- proc.time()[3]

toc - tic

