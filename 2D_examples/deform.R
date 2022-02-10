#######################################
# Visualizing the deformation: X to W #
# Steve Walsh                         #
# February 2 2022                     #
#######################################

library(Morpho)

niters <- 50000
k <- 20   # show every kth deformation
pmx <- T  # use MCMCs with Prior Mean X for hidden layer?

pdf("pdf/deform/1to18_every20th.pdf")

for (ste in 1:18) {
  
  # load in training points, deformed points, param chains, etc
  load(paste0("rda/FL_fits/storm",ste,"_niters",niters,
              "krig",if(pmx){"pmx"},".rda"))
  
  # the deformation often negates some of the values; remove this
  fit$ww <- list()
  for (i in 1:length(fit$g)) {
    fit$ww[[i]] <- ifelse(fit$w[[i]] < 0, -fit$w[[i]], fit$w[[i]])
  }
  
  # # Basic comparison of the lat/lon and deformed locations 
  # pdf("pdf/deform.pdf")
  # par(mfrow=c(2,2))
  # plot(fit$x)
  # plot(rasterFromXYZ(cbind(fit$x, fit$y)))
  # plot(fit$ww[[5]])
  # plot(fit$ww[[50000]])
  # 
  # par(mfrow=c(1,1))
  # for (i in 1:(100)) { plot(fit$ww[[i*500]], main = i*500) }
  # dev.off()
  
  
  # These two show a grid warped based on thin plate splines
  
  # # focus on deformations early on in the chain; every k'th
  # pdf(paste0("pdf/deform/",ste,if(pmx){"pmx"},".pdf"))
  for(i in 1:20){ deformGrid2d(fit$x, fit$ww[[k*(i-1)+1]],
                               ngrid=25, pch=19, main=paste(ste, k*(i-1)+1), gridcol = "black") }
  # dev.off()
  
  # look at more of big picture changes in the chain; 20 images over the chain
  # pdf(paste0("pdf/deform/",ste,if(pmx){"pmx"},"_all.pdf"))
  # for(i in 1:20){ deformGrid2d(fit$x, fit$ww[[niters/20*(i-1)+1]],
  #                              ngrid=25, pch=19, main=paste(ste, niters/20*(i-1)+1), gridcol = "black") }
  # dev.off()

}

dev.off()
