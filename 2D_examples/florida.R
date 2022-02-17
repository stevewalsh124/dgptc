######################################
# Modify FL storms: Proof of Concept #
# Steve Walsh                        #
# February 3 2022                    #
######################################

library(raster)
library(stringr) #str_locate_all 

# We have 47 training storms, 21 of which were labeled as FL landfall.
# Of these 21 TCs, 3 storms (5,10,29) made most precip elsewhere
# So, work with the remaining 18 storms, and mask out grid points
# far from the Florida peninsula for common lat/lon & proof of concept.

FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]

tc_files <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/",
                       full.names = T)
# Do this once:
# Use storm 3's points as a mask to ignore more inland grid points
# tc <- read.csv(tc_files[3], row.names = 1)
# plot(rasterFromXYZ(tc[,c(2,3,1)]), main="storm 3, to be mask")
# mask_FL <- rasterFromXYZ(tc[,c(2,3,1)])
# values(mask_FL) <- ifelse(is.na(values(mask_FL)), NA, 1)
# plot(mask_FL)
# save(mask_FL, file = "rda/FL_mask")
load("rda/FL_mask")

# ########################################################################
# # Do this once to make the FL_df (saved as rda) files in FL_storms dir #
# # Each buffered with FL_mask above in order to have same n per TC      #
##########################################################################
# storm.dirs <- list.dirs("~/NAMandST4", recursive = F)
# storms.out.of.hurdat <- namevec <- c()
# NAMlist <- ST4list <- list()
# for(i in fl){
#         # i <- 1
#         currentST4_NAM_folders <- list.dirs(storm.dirs[i])[-1] #[-1] avoids parent folder
#         ST4_folder <- currentST4_NAM_folders[1]
#         NAM_folder <- currentST4_NAM_folders[2] #uncomment the [2] as well
#         storm_yearname <- list.dirs(storm.dirs[i], full.names = F)[-1][1]
#         
#         name_start <- tail(str_locate_all(pattern ='/', ST4_folder)[[1]],1)
#         storm_name <- substr(ST4_folder,name_start[1]+5,nchar(ST4_folder))
#         namevec[i] <- storm_name
#         print(paste0("Storm #",i,": ",storm_name))
#         
#         hurdat_csv <- list.files(storm.dirs[i], pattern = ".csv")
#         
#         print(length(list.files(NAM_folder))) #should be 15 or 18, depending if 72hr timestep exists
#         NAMlist[[i]] <- list()
#         ST4list[[i]] <- list()
#         
#         # figure out if hour started is 00/12, or 06/18
#         hr6 <- as.numeric(substr(list.files(ST4_folder)[4],13,14)) %% 12 == 6
#         # if(!hr6) next
#         #Identifying where in the NAM folder the first 24 hours after landfall files are (latest timestep)
#         if(length(list.files(NAM_folder, full.names = T))%in%c(21,29) & hr6){first24 <- 2:9} else {
#                 if(length(list.files(NAM_folder, full.names = T))==53 & hr6){first24 <- c(4,7,10,13,16,19,22,25)} else { #2:25
#                         if(length(list.files(NAM_folder, full.names = T))==21 & !hr6){first24 <- c(5,9)} else {
#                                 if(length(list.files(NAM_folder,full.names = T))==18){first24 <- c(7,8)} else {
#                                         if(length(list.files(NAM_folder,full.names = T))==15){first24 <- c(5,6)} else {
#                                                 if(length(list.files(NAM_folder,full.names = T))==12&storm_name=="rita"){first24 <- c(5,6)} else {
#                                                         if(length(list.files(NAM_folder, full.names = T))==53){first24 <- c(13,25)} else {
#                                                                 print("Total NAM files for this storm not 15 or 18 or 53; skipping"); next}}}}}}}
#         
#         for (k in 1:length(list.files(ST4_folder,full.names = T)[1:4])) {
#                 # k <-1
#                 ST4file <- list.files(ST4_folder,full.names = T)[k]
#                 ST4list[[i]][[k]] <- (raster(ST4file) %>% projectRaster(crs = "+proj=longlat +datum=WGS84", method = "ngb"))
#                 ST4list[[i]][[k]] <- raster::resample(ST4list[[i]][[k]], mask_FL, method = "ngb")
#                 
#                 print(min(values(ST4list[[i]][[k]]), na.rm=T))
#                 # values(ST4list[[i]][[k]])[which(values(ST4list[[i]][[k]]<=1))] <- 1 
#                 #change nonpositive projected precipitation to 0.0000000001 to implement log precip wihtout -Inf
#         }   
#         
#         
#         for (j in 1:length(list.files(NAM_folder,full.names = T)[first24])) {
#                 # j <-1
#                 NAMfile <- list.files(NAM_folder,full.names = T)[first24[j]]
#                 
#                 #Establish correct NAM band for total precip, mostly depends on date of storm (3 w issues)
#                 date_start <- str_locate_all(pattern ='_218_', NAMfile)[[1]][1,2] + 1
#                 storm_date <- as.numeric(substr(NAMfile,date_start,date_start+7))
#                 if(storm_date < 20060701){NAM_band <- 10} else 
#                         if(storm_date > 20060701 && storm_date < 20170101){NAM_band <- 11} else {NAM_band <- 373}
#                 print(NAM_band)
#                 NAMlist[[i]][[j]] <- (raster(NAMfile, band=NAM_band) %>% projectRaster(crs = "+proj=longlat +datum=WGS84",
#                                                                                        method = "ngb"))
#                 print(min(values(NAMlist[[i]][[j]]), na.rm=T))
#                 # values(NAMlist[[i]][[j]])[which(values(NAMlist[[i]][[j]]<=1))] <- 1 
#                 #change nonpositive projected precipitation to 0.0000000001 to implement log precip without -Inf
#         }
#         
#         ST4_first12  <- Reduce("+",ST4list[[i]])*mask_FL
#         ST4_df_first12 <- rasterToPoints(ST4_first12)
#         
#         NAM_first12 <- Reduce("+",NAMlist[[i]])*mask_FL
#         NAM_df_first12  <- rasterToPoints(NAM_first12)
#         
#         NAM_plotter <- sqrt(NAM_first12)
#         ST4_plotter <- sqrt(ST4_first12)
#         
#         ## Find the extra points from the files and remove them
#         ## When NAM and ST4 have different amounts of pixels
#         ## if NAM has more rows...
#         if(nrow(rasterToPoints(NAM_plotter)) > nrow(rasterToPoints(ST4_plotter))){
#                 ST40 <- ST4_plotter
#                 values(ST40)[!is.na(values(ST40))] <- 0
#                 NAM_plotter <- NAM_plotter - ST40
#         }
#         
#         ## if ST4 has more rows
#         if(nrow(rasterToPoints(NAM_plotter)) < nrow(rasterToPoints(ST4_plotter))){
#                 NAM0 <- NAM_plotter
#                 values(NAM0)[!is.na(values(NAM0))] <- 0
#                 ST4_plotter <- ST4_plotter - NAM0
#         }
#         
#         if(nrow(rasterToPoints(NAM_plotter)) > nrow(rasterToPoints(ST4_plotter))){
#                 ST40 <- ST4_plotter
#                 values(ST40)[!is.na(values(ST40))] <- 0
#                 NAM_plotter <- NAM_plotter - ST40
#         }
#         
#         error <- ST4_plotter - NAM_plotter
#         FL_df <- rasterToPoints(error)
#         save(FL_df, file = paste0("rda/FL_storms/",i,".rda"))
#         
# }
##########################################################################

# Plot each FL TC and check for same number of grid points
pdf("pdf/FL_storms.pdf")
par(mfrow=c(2,3))
for (ste in fl) {
        load(paste0("rda/FL_storms/",ste,".rda")) 
        print(c(ste, nrow(FL_df)))
        plot(rasterFromXYZ(FL_df)*mask_FL, main=ste)
}
dev.off()
