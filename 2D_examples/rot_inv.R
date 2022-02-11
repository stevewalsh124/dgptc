#############################################
# Make rotation-invariant deformation plots #
#############################################

pdf("pdf/deform/deform_compare.pdf")

for (ste in c(1,3,5,6,7,8,10,11,12,13,15,16,17)) {
  load(paste0("~/dgptc/2D_examples/rda/FL_fits/storm",ste,"_niters5000krigpmx.rda"))
  
  X <- fit$x
  
  k <- 100 # every kth warping to plot
  tp <- 10 # total plots to produce 
  
  # the deformation often negates some of the values; remove negatives
  fit$ww <- list()
  for (i in 1:length(fit$g)) {
    fit$ww[[i]] <- ifelse(fit$w[[i]] < 0, -fit$w[[i]], fit$w[[i]])
  }
  W <- fit$ww
  
  library(Morpho)
  
  # The original non-invariant plots; shown in left side below for comparison
  for (i in 1) { #1:tp
    deformGrid2d(fit$x, fit$ww[[k*(i-1)+1]],
                 ngrid=25, pch=19, main=paste(ste, k*(i-1)+1), gridcol = "black")
  }
  
  # The Xs beforehand
  # plot(X, asp = 1, pch=".")
  # plot(W[[4]], asp = 1, pch=".")
  
  # (1) Call a point in the southern border of Florida point (0,0);
  # in this storm, find first (and only for this TC) point with y=0; this is (0.8,0)
  # So I will subtract 0.8 from each of these
  lo_lat <- which(X[,2]==0 & abs(X[,1]-0.7692308) < 1e-3)#which(X[,2]==0)[1]
  X[lo_lat,]
  
  # (2) Find the point on the northern border of Florida with same longitude as (0,0) and call it point (0,1);
  hi_lat <- which(X[,2]==1 & abs(X[,1]-0.7692308) < 1e-3)
  X[hi_lat,]
  
  # (c) Find the translation that makes (x_1,x_2) into (0,0) (that's easy). 
  x_adj <- X[lo_lat,1]
  X[,1] <- X[,1] - x_adj
  
  # (3) For each simulated hidden layer: 
  for(i in k*(1:tp)+1){
    par(mfrow=c(1,2))
    
    deformGrid2d(fit$x, fit$ww[[i]],
                 ngrid=25, pch=19, main=paste(ste, i, "old"), gridcol = "black")
    
    # i <- 1 (above)
    #   (a) Find the "deformed point" corresponding to point (0,0), say (x_1,x_2);
    #   Note: Ws have not had x_adj subtracted from longitudes
    w <- W[[i]]
    w[lo_lat,] 
    
    # (b) Find the "deformed point" corresponding to point (0,1), say (x_3,x_4);
    w[hi_lat,]
    
    #     (d) Apply that translation to all deformed points. For example, after (x_3,x_4) 
    #           gets translated, it becomes (x_3-x_1, x_4-x_2);
    w[,1] <- w[,1] - x_adj
    # plot(w, asp=1, pch=".")
    which(X[,1]==0 & X[,2]==0)
    
    #     (e) Now use the dot product to compute the angle between the vector that 
    #           goes from (0,0) to (x_3-x_1, x_4-x_2). Find the corresponding 
    #           rotation-scale matrix that makes (x_3-x_1, x_4-x_2) into (0,1);
    u <- X[hi_lat,]
    v <- w[hi_lat,]
    
    # find angle between (0,0) and v will result in arccos(0/0)=NaN
    # find angle between ~(0,1) and its corresponding deformation 
    ang <- function(u,v){acos((t(u)%*%v)/((sqrt(t(u)%*%u) * sqrt(t(v)%*%v))))}
    
    phi <- ang(u,v)
    rotMtx <- matrix(c(cos(phi),sin(phi),-sin(phi),cos(phi)),2,2)
    
    # rotMtx%*%u / v should return two numbers of equal value
    # they don't, so I tried taking mean, or the max(abs())
    print(i)
    print(c(rotMtx%*%u / v))
    r <- abs(max(rotMtx%*%u / v))
    
    #      (f) Apply that rotation-scale matrix to all translated points. 
    f <- function(u){r*rotMtx %*%u}
    ww <- t(apply(w, 1, f))
    ww[,1] <- ww[,1]+x_adj
    # plot(ww)
    #      (g)  I think that will give us a deformation that is rotation-invariant.
    deformGrid2d(fit$x, ww,
                 ngrid=25, pch=19, main=paste(ste, i), gridcol = "black")
  }
}

dev.off()
