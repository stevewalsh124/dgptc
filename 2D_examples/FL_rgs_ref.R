###################################
# Get range of lat/lon for FL TCs #
# Also establish two reference    #
# points for deformation          #
#                                 #
# Steve Walsh                     #
# February 13 2022                #
###################################

# First, find range of lat and lon to scale
# all points to be in the unit square.
# Save these points as FL_rgs

# Second, use FL_rgs to rescale all FL TCs,
# and then use the intersection of all these
# collections of grid points to find a 
# common set of reference points
# Save these points as FL_ref

library(raster)

FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]

# # For transformation to [0,1]; saved below
# # Init mins and maxes
# min_x <- min_y <- Inf
# max_x <- max_y <- -Inf
load("rda/FL_rgs.rda")

pdf("pdf/intersection_FLs.pdf")
for (ste in fl) {
  # loads FL_df
  load(paste0("~/dgptc/2D_examples/rda/FL_storms/",ste,".rda"))

  # Read in the specific storm data frame
  tc <- data.frame(FL_df)
  names(tc) <- c("x","y","value")
  
  # scaled lon and lat to be in [0,1]
  tc$xs <- (tc$x - FL_rgs[1])/(FL_rgs[2] - FL_rgs[1])
  tc$ys <- (tc$y - FL_rgs[3])/(FL_rgs[4] - FL_rgs[3])  
  
  # # Get min and max that is uniform over all (FL) TCs
  # # For transformation to [0,1]; saved below
  # min_x <- min(min_x, FL_df[,"x"])
  # max_x <- max(max_x, FL_df[,"x"])
  # min_y <- min(min_y, FL_df[,"y"])
  # max_y <- max(max_y, FL_df[,"y"])
  
  # Plot the intersection of the all TCs up to the ste'th
  r <- raster::rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value))
  if(ste == fl[1]){r_sum <- r} else {r_sum <- r_sum + r}
  plot(r_sum, main = ste)
}
dev.off()

# # Did this once; saved results to get FL TCs to common [0,1]
# FL_rgs <- c(min_x, max_x, min_y, max_y)
# save(FL_rgs, file = "rda/FL_rgs.rda")

# Find the (0,0) zero zero point (zz)
r_sum_df <- rasterToPoints(r_sum)
wl <- which(r_sum_df[,"y"] == min(r_sum_df[,"y"]))
zz <- r_sum_df[wl,1:2]

# only look at x (aka lon) value that has the (0,0) point 
# points order top to bottom, so (0,1) zero-one (zo) point 
# is 1st in the subset
r_sub <- r_sum_df[which(r_sum_df[,"x"]==r_sum_df[wl,"x"]),]
zo <- r_sub[1,1:2]
FL_ref <- rbind(zz,zo)
save(FL_ref, file = "rda/FL_ref.rda")
