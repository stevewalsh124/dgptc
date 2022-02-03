#######################################
# Visualizing the deformation: X to W #
# Steve Walsh                         #
# February 2 2022                     #
#######################################

# load in training points, deformed points, param chains, etc
load("rda/storm11_niters10000.rda")

# the deformation often negates some of the values; remove this
fit$ww <- list()
for (i in 1:length(fit$g)) {
  fit$ww[[i]] <- ifelse(fit$w[[i]] < 0, -fit$w[[i]], fit$w[[i]])
}

# Basic comparison of the lat/lon and deformed locations 
pdf("pdf/deform.pdf")
par(mfrow=c(2,2))
plot(fit$x)
plot(rasterFromXYZ(cbind(fit$x, fit$y)))
plot(fit$ww[[5]])
plot(fit$ww[[50000]])

par(mfrow=c(1,1))
for (i in 1:(100)) { plot(fit$ww[[i*500]], main = i*500) }

dev.off()


# These two show a grid warped based on thin plate splines

# focus on deformations early on in the chain; every 50th
pdf("pdf/deform_grid.pdf")
for(i in 1:20){ deformGrid2d(fit$x, fit$ww[[50*(i-1)+1]],
                             ngrid=25, pch=19, main=50*(i-1)+1, gridcol = "black") }
dev.off()

# look at more of big picture changes in the chain; every 500th
pdf("pdf/deform_grid_all.pdf")
for(i in 1:20){ deformGrid2d(fit$x, fit$ww[[2500*(i-1)+1]],
                             ngrid=25, pch=19, main=2500*(i-1)+1, gridcol = "black") }
dev.off()

