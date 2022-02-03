######################################
# Modify FL storms: Proof of Concept #
# Steve Walsh                        #
# February 3 2022                    #
######################################

library(raster)

# We have 47 training storms, 21 of which were labeled as FL landfall.
# Of these 21 TCs, 3 storms (5,10,29) made most precip elsewhere
# So, work with the remaining 18 storms, and mask out grid points
# far from the Florida peninsula for common lat/lon & proof of concept.

FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]

tc_files <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/",
                       full.names = T)

# Use storm 3's points as a mask to ignore more inland grid points
# tc <- read.csv(tc_files[3], row.names = 1)
# plot(rasterFromXYZ(tc[,c(2,3,1)]), main="storm 3, to be mask")
# mask_FL <- rasterFromXYZ(tc[,c(2,3,1)])
# values(mask_FL) <- ifelse(is.na(values(mask_FL)), NA, 1)
# plot(mask_FL)
# save(mask_FL, file = "rda/FL_mask")
load("rda/FL_mask")

par(mfrow=c(2,3))
for (ste in fl) {
 tc <- read.csv(tc_files[ste], row.names = 1)
 plot(rasterFromXYZ(tc[,c(2,3,1)])*mask_FL, main=ste)
 FL_df <- rasterToPoints(rasterFromXYZ(tc[,c(2,3,1)])*mask_FL)
 save(FL_df, file=paste0("rda/FL_storms/",ste,".rda"))
}
