#################################
# Parameter estimation for TCs  #
# January 31, 2022              #
# Steve Walsh, Annie Sauer      #
#################################

library(raster)
library(deepgp) # version >= 0.3.0

tic <- proc.time()[3]

# Do you want to...
# use the FL subset TCs?
# specify no warping a priori? (ie, should W have prior mean X?)
# do any spatial prediction (kriging)?
do_FL <- T
pmx <- T
krig <- T

# Read in the previous burned-in values for params and w
# loads init_param and init_w
# load(paste0("rda/burn_params_FL",if(pmx){"pmx"},".rda"))
# load(paste0("rda/burn_w_FL",if(pmx){"pmx"},".rda"))

# Load the appropriate files (most FL storms, or all full storms)
if(do_FL){
  tc_files <- list.files("rda/FL_storms", full.names = T)
} else {
  tc_files <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/",
                         full.names = T)
}

# storm to evaluate (11 is smallest)
ste <- 11

# number of iterations for MCMC
niters <- 25001

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# Read in the specific storm data frame
if(do_FL){
  load(tc_files[ste])
  tc <- data.frame(FL_df)
  names(tc) <- c("x","y","value")
} else {
  tc <- read.csv(tc_files[ste], row.names = 1)
}

# scaled lon and lat to be in [0,1]
tc$xs <- (tc$x - min(tc$x))/(max(tc$x)-min(tc$x))
tc$ys <- (tc$y - min(tc$y))/(max(tc$y)-min(tc$y))
tc$zs <- (tc$value - mean(tc$value))/(sd(tc$value))

# set training and test sets
set.seed(1) # keep a consisten set of train/test
wl <- which(tc$ys==0 & abs(tc$xs-0.7692308) < 1e-3) # save a (0,0) point
wh <- which(tc$ys==1 & abs(tc$xs-0.7692308) < 1e-3) # save a (0,1) point
all_but_ws <- (1:nrow(tc))[-c(wl,wh)]
train <- c(sample(all_but_ws, max(500, floor(.25*nrow(tc)))),wl,wh)
test <- (1:nrow(tc))[-train]

# subset training data
tc_samp <- tc[train,]
x <- cbind(tc_samp$xs, tc_samp$ys)
y <- tc_samp$zs
# plot(rasterFromXYZ(data.frame(cbind(x,y))))

# Locations for predictions
tc_pred <- tc[test,]
xx <- cbind(tc_pred$xs, tc_pred$ys)

# Get appropriate HURDAT info depending on FL or not
FL <- read.csv("~/NAM-Model-Validation/csv/stormsFL.csv", row.names = 1)$x
fl <- FL[-which(FL %in% c(5,10,29))]
if(do_FL){ind <- fl[ste]} else {ind <- ste}

# Read in each storm's hurdat info (track information)
# 6-hourly info that was interpolated to have 3-hourly info (by Stephanie Zick)
new.spline <- T
if(new.spline){
  hurdat <- read.csv(list.files( "/home/walsh124/NAM-Model-Validation/csv/BestTrackSpline15min/",
                                 pattern = ".csv", full.names = T)[ind], row.names = 1)
  colnames(hurdat) <- c("DateTime", "Hour", "Lon", "Lat", "MSLP", "Winds..knots.")
  if(hurdat$Lon[1] > 0) hurdat$Lon <- hurdat$Lon-360
} else {
  hurdat <- read.csv(list.files(storm.dirs[ind], pattern = ".csv", full.names = T))
  colnames(hurdat) <- c("BTid", "SNum", "Year", "Month", "Day", "Hour", "Lat", "Lon", "Winds", "Pressure")
}

## Get the subset of the hurdat track info corresponding to the 24 hour period
# ST4 folder is the 2nd listed in the storm dir
storm.dirs <- list.dirs("/home/walsh124/NAMandST4", recursive = F)
ST4_folder <- list.dirs(storm.dirs[ind])[2]
eye_1 <- list.files(ST4_folder,full.names = T)[1]
library(stringr)
time_start  <- tail(str_locate_all(pattern ='/', eye_1)[[1]],1)[,2]
storm_time  <- substr(eye_1,time_start[1]+5,nchar(eye_1))
storm_year  <- substr(storm_time,1,4)
storm_month <- substr(storm_time,5,6)
storm_day   <- substr(storm_time,7,8)
storm_hour  <- substr(storm_time,9,10)

print(paste("Storm # is", ste, "Year is",storm_year,"Month is",storm_month,
            "Day is",storm_day,"Hour is",storm_hour))

raster_locs <- my_NS_nsv$coords
if(kern.locs) raster_locs <- my_NS_nsv$mc.locations

if(new.spline){
  ind <- which(substr(hurdat$DateTime, 9, 13)==paste(storm_day, storm_hour))[1]
  hurdat_sub <- hurdat[(ind-24):(ind+72), ] # -12, +36 for 30 min
  hurdat_sub <- hurdat_sub[!is.na(hurdat_sub$Hour),]
  # Get the hurricane track locations and the raster locations
  track_locs <- hurdat_sub[,3:4]
} else {
  ind <- which(hurdat$Day==as.integer(storm_day) & hurdat$Hour==as.integer(storm_hour))
  hurdat_sub <- hurdat[(ind-2):(ind+5), ]#7:8]
  hurdat_sub <- hurdat_sub[!is.na(hurdat_sub$Hour),]
  # Get the hurricane track locations and the raster locations
  track_locs <- hurdat_sub[,8:7]
}

# Calculate the SQUARED distances between the track and raster locations
# Obtain the minimum distance from the track for each raster grid point
if(orig_dist){
  dists <- sqrt(plgp::distance(raster_locs, track_locs))
  coast_dists <- sqrt(plgp::distance(raster_locs, sea_locs))
} else {
  dists <- plgp::distance(raster_locs, track_locs)
  coast_dists <- plgp::distance(raster_locs, sea_locs)
}

min_dists <- apply(dists, 1, min)
min_coast_dists <- apply(coast_dists, 1, min)


# Fit two-layer DGP (exponential cov fn)
if(pmx){
  fit <- fit_two_layer(x, y, nmcmc = niters, cov = "matern", v=0.5, vecchia = T, 
                       # theta_y_0 = init_param[ste,1], 
                       # theta_w_0 = init_param[ste,2:3], 
                       # w_0 = init_w[[ste]][[1]],
                       true_g = sqrt(.Machine$double.eps),
                       settings = list(w_prior_mean = x))
} else {
  fit <- fit_two_layer(x, y, nmcmc = niters, cov = "matern", v=0.5, vecchia = T, 
                       # theta_y_0 = init_param[ste,1], 
                       # theta_w_0 = init_param[ste,2:3], 
                       # w_0 = init_w[[ste]][[1]],
                       true_g = sqrt(.Machine$double.eps))
}

# fit <- trim(fit, 1000, 1) # retain 2500 samples

# save before predict
if(do_FL){
  save(fit, file = paste0("rda/FL_fits/storm",ste,"_niters",niters,".rda"))
} else {
  save(fit, file = paste0("rda/storm",ste,"_niters",niters,".rda"))
}

# predict and save after
if(krig){
  fit <- predict(fit, xx)
  save(fit, file = paste0("rda/",if(do_FL){"FL_fits/"},
                          "storm",ste,
                          "_niters",niters,
                          if(krig){"krig"},
                          if(pmx){"pmx"},".rda"))
}

if(krig){
  # Combine results
  pred <- data.frame(xx = xx, mean = fit$mean, s2 = fit$s2_smooth)
  
  # Save prediction results
  write.csv(pred, paste0("csv/",if(do_FL){"FL_fits/"},"storm",ste,"_niters",
                         niters, if(pmx){"pmx"},".csv"), row.names = FALSE)
  
  # get range for plots
  rg <- range(c(tc$value, pred$mean)) #, pred2$mean
  s2_rg <- range(c(pred$s2)) #, pred2$s2
}

# compare predictions via plots
pdf(paste0("pdf/",if(do_FL){"FL_fits/"},"storm",ste,
           "_niters_",niters,if(krig){"krig"},if(pmx){"pmx"},".pdf"))
if(krig){par(mfrow=c(1,2))}
plot(rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value)), if(krig){zlim=rg}, main="EF")
if(krig){
  plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$mean)), zlim=rg, main="pred mean")
  # plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$mean)), zlim=rg)
  
  plot(rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value)), main="EF")
  plot(rasterFromXYZ(cbind(pred$xx.1,pred$xx.2,pred$s2)), zlim=s2_rg, main="pred s2")
  # plot(rasterFromXYZ(cbind(pred2$xx.1,pred2$xx.2,pred2$s2)), zlim=s2_rg)
}

par(mfrow=c(1,3))
plot(fit$theta_w[,1], type="l", main="theta_w[,1]")
plot(fit$theta_w[,2], type="l", main="theta_w[,2]")
plot(fit$theta_y, type="l", main="theta_y")

plot(fit$tau2, type="l", main="tau2")
plot(fit$g, type="l", main="g")

# plot(fit2$theta_w[,1], type="l")
# plot(fit2$theta_w[,2], type="l")
# plot(fit2$theta_y, type="l")
dev.off()

toc <- proc.time()[3]

toc - tic
