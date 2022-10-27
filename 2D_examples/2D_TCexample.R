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
# do any spatial prediction (kriging)? Note: this is not needed for UQ.
do_FL <- T
pmx <- F
krig <- F
vecchia <- T

cov <- "matern"
v <- 1.5
true_g = NULL #sqrt(.Machine$double.eps)

# what to evaluate
if(do_FL){
  wte <- "error" # "NAM" "ST4"
  if(!(wte %in% c("error","NAM","ST4"))) stop("wte should be %in% c('error','NAM','ST4')")
}

# Read in the previous burned-in values for params and w
# loads init_param and init_w
# load(paste0("rda/burn_params_FL",if(pmx){"pmx"},".rda"))
# load(paste0("rda/burn_w_FL",if(pmx){"pmx"},".rda"))

# Load the appropriate files (most FL storms, or all full storms)
if(do_FL){
  tc_files <- list.files(paste0("rda/FL_storms",if(wte=="NAM"){"/NAM"},if(wte=="ST4"){"/ST4"}), 
                         full.names = T, pattern = ".rda")
} else {
  tc_files <- list.files("~/NAM-Model-Validation/csv/error_df/subtractPWmeanF/",
                         full.names = T)
}

# storm to evaluate (11 is smallest if FL_df <- FALSE)
ste <- 18

# number of iterations for MCMC
niters <- 12000
nburn <- 2000
kth <- 5

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

# Read in the specific storm data frame
if(do_FL){
  load(tc_files[ste])
  if(wte == "error") tc <- data.frame(FL_df)
  if(wte == "NAM"){
    NAM_df <- rasterToPoints(NAM_plotter)
    tc <- data.frame(NAM_df)
  }
  if(wte == "ST4"){
    ST4_df <- rasterToPoints(ST4_plotter)
    tc <- data.frame(ST4_df)
  }
  names(tc) <- c("x","y","value")
} else {
  tc <- read.csv(tc_files[ste], row.names = 1)
}

# scaled lon and lat to be in [0,1]
load("rda/FL_rgs.rda")
tc$xs <- (tc$x - FL_rgs[1])/(FL_rgs[2] - FL_rgs[1])
tc$ys <- (tc$y - FL_rgs[3])/(FL_rgs[4] - FL_rgs[3])
tc$zs <- (tc$value - mean(tc$value))/(sd(tc$value))

# set training and test sets
set.seed(1) # keep a consisten set of train/test
load("rda/FL_ref.rda")
zz <- FL_ref[1,]
zo <- FL_ref[2,]
sw <- FL_ref[3,]
se <- FL_ref[4,]
wl <- which(abs(tc$xs-zz[1]) < 1e-3 & abs(tc$ys-zz[2]) < 1e-3) # (0,0) point
wh <- which(abs(tc$xs-zo[1]) < 1e-3 & abs(tc$ys-zo[2]) < 1e-3) # (0,1) point
ww <- which(abs(tc$xs-sw[1]) < 1e-3 & abs(tc$ys-sw[2]) < 1e-3) # (0,1) point
we <- which(abs(tc$xs-se[1]) < 1e-3 & abs(tc$ys-se[2]) < 1e-3) # (0,1) point
all_but_ws <- (1:nrow(tc))[-c(wl,wh,ww,we)]
train <- c(sample(all_but_ws, 496),wl,wh,ww,we)
test <- (1:nrow(tc))[-train]

# subset training data
tc_samp <- tc[train,]
x <- cbind(tc_samp$xs, tc_samp$ys)
y <- tc_samp$zs
# plot(rasterFromXYZ(data.frame(cbind(x,y))))

# Locations for predictions
tc_pred <- tc[test,]
xx <- cbind(tc_pred$xs, tc_pred$ys)

# Fit two-layer DGP (exponential cov fn)
if(pmx){
  fit <- fit_two_layer(x, y, nmcmc = niters, cov = cov, v = v, vecchia = vecchia, 
                       # theta_y_0 = init_param[ste,1],
                       # theta_w_0 = 1000/(.0001/1000),
                       # w_0 = init_w[[ste]][[1]],
                       true_g = true_g,
                       settings = list(w_prior_mean = x))#,
                                       # alpha = list(g = 1.5, theta_w = 1000, theta_y = 1.5),
                                       # beta = list(g = 3.9, theta_w = .0001/1000, theta_y = 3.9/6)))
} else {
  fit <- fit_two_layer(x, y, nmcmc = niters, cov = cov, v = v, vecchia = vecchia, 
                       # theta_y_0 = init_param[ste,1],
                       # theta_w_0 = init_param[ste,2:3],
                       # w_0 = init_w[[ste]][[1]],
                       true_g = true_g)
}

fit <- deepgp::trim(fit, nburn, kth) # retain 2500 samples

# save before predict
if(do_FL){
  save(fit, file = paste0("rda/FL_fits/storm",ste,"_niters",niters,if(pmx){"pmx"},"_",wte,".rda"))
} else {
  save(fit, file = paste0("rda/storm",ste,"_niters",niters,if(pmx){"pmx"},".rda"))
}

# predict and save after
if(krig){
  fit <- predict(fit, xx)
  save(fit, file = paste0("rda/",if(do_FL){"FL_fits/"},
                          "storm",ste,
                          "_niters",niters,
                          "krig",if(pmx){"pmx"},
                          "_",wte,".rda"))
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
           "_niters_",niters,if(krig){"krig"},if(pmx){"pmx"},"_",wte,".pdf"))
if(krig){par(mfrow=c(1,2))}
plot(rasterFromXYZ(cbind(tc$xs, tc$ys, tc$value)), if(krig){zlim=rg}, main=paste(wte, niters, cov, v, true_g))
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
