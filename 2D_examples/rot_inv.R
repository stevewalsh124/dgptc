#############################################
# Make rotation-invariant deformation plots #
# Steve Walsh                               #
# February 14 2022                          #
#############################################

library(Morpho) #deformGrid2d
library(marmap) #griddify

k <- 24999 # every kth warping to plot
tp <- 1 # total plots to produce (excluding original)

FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]

pdf(paste0("pdf/deform/deform_compare_k",k,"_tp",tp,".pdf"))

load("rda/all_FL_tracks.rda")
load("rda/FL_rgs.rda")
for (i in fl) {
  all_FL_tracks[[i]][,1] <- (all_FL_tracks[[i]][,1] - FL_rgs[1]) / 
    (FL_rgs[2] - FL_rgs[1])
  all_FL_tracks[[i]][,2] <- (all_FL_tracks[[i]][,2] - FL_rgs[3]) / 
    (FL_rgs[4] - FL_rgs[3])
}

# Loop over TCs
for (ste in 1:18) {
  load(paste0("rda/FL_fits/storm",ste,"_niters25000.rda")) #loads fit
  
  x <- t(fit$x)
  
  # the deformation often negates some of the values; remove negatives
  # this removes any reflections incurred within the hidden layer
  W <- list()
  for (i in 1:fit$nmcmc) W[[i]] <- ifelse(fit$w[[i]] < 0, -fit$w[[i]], fit$w[[i]])

  # The original non-invariant plots; shown in left side for comparison
  # par(mfrow=c(1,1))
  # deformGrid2d(fit$x, W[[1]], ngrid=25, pch=19, main=paste(ste, 1), gridcol = "black")
  
  # Get reference points and scale
  load("rda/FL_ref.rda")
  zz <- FL_ref[1,]
  zo <- FL_ref[2,]
  # points(zz[1],zz[2], lwd=3) # add to plot above
  # points(zo[1],zo[2], lwd=3) # add to plot above
  
  # this finds the length of the original reference vector
  zzo <- zz - zo
  ref.scale <- sqrt(t(zzo)%*%zzo)
  
  # establish which points in x (and y) are reference points
  # which low (wl) is (0,0) and which high (wh) is (0,1)
  wl <- which(abs(x[1,]-FL_ref[1,1]) < 1e-3 & abs(x[2,]-FL_ref[1,2]) < 1e-3) # (0,0) point
  wh <- which(abs(x[1,]-FL_ref[2,1]) < 1e-3 & abs(x[2,]-FL_ref[2,2]) < 1e-3) # (0,1) point
  
  # If reference point pair are not in training set, skip
  if(length(wl)==0 | length(wh)==0) next
  print(ste)
  
  # For each simulated hidden layer: 
  for(i in k*(1:tp)+1){
    
    par(mfrow=c(1,2))
    # On the left: the original (unchanged) deformation
    deformGrid2d(fit$x, W[[i]], ngrid=25, pch=19, main=paste(ste, i, "old"), gridcol = "black")
    
    # Take these deformed points (y) and rotate/scale/translate them
    # to match orig lat/lon points (x) as well as possible
    y <- t(W[[i]])

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
    ww <- t(y.rot.scale.tran)
    
    # On right: plot a deformation that is rotation-invariant wrt original lat/lon
    deformGrid2d(fit$x, ww, ngrid=25, pch=19, main=paste(ste, i), gridcol = "black")
    
    yy <- fit$y
    irreg <- as.data.frame(cbind(ww,yy))
    colnames(irreg) <- c("lon","lat","y")
    # use griddify to create a 40x60 grid
    reg <- griddify(irreg, nlon = 40, nlat = 60)
    
    # Plot the new bathy object and overlay the original data points
    plot(fit$x[,1], fit$x[,2], pch = ".", cex = 0.3, asp=1,
         xlim = c(range(fit$x[,1], irreg$lon)), ylim = c(range(fit$x[,2], irreg$lat)))
    plot(reg, add=T, alpha=1)
    points(irreg$lon, irreg$lat, pch = ".", cex = 0.3, col = col2alpha(1, alpha = 0.1))
    
    
    plot(rasterFromXYZ( cbind(fit$x, fit$y)))
    points(all_FL_tracks[[fl[ste]]][,1], all_FL_tracks[[fl[ste]]][,2], pch=".")
    
  }
}

dev.off()
