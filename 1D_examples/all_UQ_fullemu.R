# coompare UQ results for each of 111 model batches
# with/without mean subtracted out
# also compare 111 warpings

zero_mean <- T
do_warps <- F
mtes <- 111

PDF <- T
if(PDF) pdf(paste0("pdf/all_UQ_fullemu",if(zero_mean){"_zm"},".pdf"))

vecchia <- F
pmx <- F
one_layer <- F
force_id_warp <- F

# Smooth the precision info so it's not a step function for each data type (pt, lo, hi)?
smooth_precs <- T
if(smooth_precs) k_sm <- 10 #rolling mean uses k numbers

# Taper the covariance matrix before the MCMC fit?
taper_cov <- F

# Do a kriging step?
krig <- F

# Use hi res in Ybar calculation?
use_hi <- T

# Model the correlated errors with a covariance function?
cf_errors <- T
if(cf_errors){
  err_cov <- "matern"#"exp2"#
  loess_span <- 0.15
  err_v   <- paste0("est", loess_span)#ifelse(err_cov == "matern", 2.5, 999)
  err_g   <- NULL#sqrt(.Machine$double.eps)#
  err_g_msg <- ifelse(is.null(err_g),"estg","fixg")
}

ncores <- 2
tolpower <- -10

cov_fn <- "matern"#"exp2"#

if(taper_cov) tau_b <- .2
nrun <- 16
nmcmc <- 5000
nburn <- 2000
kth <- 4

if(zero_mean){
  for (mte in 1:mtes){
      print(mte)
    load(paste0("/projects/precipit/1D_real_study/rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))
    
    m <- fitcov$m
    
    # png(paste0("png/model",mte,"_emuspace_rmAvg.png"), width=4000, height = 2400, res=400)
    par(mar=c(4,4.5,1,1), mfrow=c(1,1))
    plot(fitcov$x, fitcov$y - m, type="n", #xlim = log10(c(.04,.35)),
         ylim = c(-.33,.33),
         xlab=expression(log[10](k)),
         ylab='script P', main = mte)#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
    for (i in 1:16) lines(fitcov$x, Y[i,] - m, col="gray", lwd=1)
    lines(fitcov$x, colMeans(Y) - m, col=cb_cols[2], lwd=2)
    lines(fitcov$x, y_hi - m, col=cb_cols[3], lwd=2)
    # lines(fitcov$x, fitcov$Ms[,1], col="red", lwd=2)
    lines(fitcov$x, fitcov$ubb - m, col=cb_cols[4], lwd=2, lty=3)
    lines(fitcov$x, fitcov$lbb - m, col=cb_cols[4], lwd=2, lty=3)
    lines(fitcov$x, cosmscrPsz - m, col=cb_cols[8], lwd=2, lty=2)
    legend("topleft",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
           col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
  }
} else {
  for (mte in 1:mtes) {
    print(mte)
    load(paste0("/projects/precipit/1D_real_study/rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))
    
    par(mar=c(4,4.5,1,1), mfrow=c(1,1))
    plot(fitcov$x, fitcov$y, type="n", #xlim = log10(c(.04,.35)),
         ylim = c(-3,2),#range(Y),
         xlab=expression(log[10](k)),
         ylab='script P', main = paste("Model batch",mte))#TeX(r'($log_{10}(k^{1.5}P(k)/2\pi^2)$)'))
    for (i in 1:16) lines(fitcov$x, Y[i,], col="gray", lwd=1)
    lines(fitcov$x, colMeans(Y), col=cb_cols[2], lwd=2)
    lines(fitcov$x, y_hi, col=cb_cols[3], lwd=2)
    # lines(fitcov$x, fitcov$Ms[,1], col="red", lwd=2)
    lines(fitcov$x, fitcov$ubb, col=cb_cols[4], lwd=2, lty=3)
    lines(fitcov$x, fitcov$lbb, col=cb_cols[4], lwd=2, lty=3)
    lines(fitcov$x, cosmscrPsz, col=cb_cols[8], lwd=2, lty=2)
    legend("bottomright",legend = c("low res","low res avg","hi res", "UQ","cosmicEMU"), 
           col = c("gray",cb_cols[c(2,3,4,8)]), lty=c(1,1,1,3,2), lwd=2)
    
  }
}

if(PDF) dev.off()

if(do_warps & PDF){
  pdf("pdf/all_Ws_fullemu.pdf")
  for (mte in 1:mtes) {
    print(mte)
    load(paste0("/projects/precipit/1D_real_study/rda/emuspace_",nmcmc,"_",one_layer,if(pmx){"_pmx"},
                if(cf_errors){paste0("_cfe",err_v,err_g_msg)},if(taper_cov){paste0("tpr",tau_b)},
                if(force_id_warp){"_fiw"},if(vecchia){"_vec"},"model",mte,"_",k_sm,".rda"))
    plot.warp(fitcov)
    mtext(mte) 
  }
  dev.off()
}
