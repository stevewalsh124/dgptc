################################
# Compare chains b/w FL storms #
# Steve Walsh                  #
# Febraury 4 2022              #
################################

pdf("pdf/all_FL_chains_25to50k.pdf")

# Choose based on the files you want to load
niters <- 50000

pmx_means <- pmx_medians <- pm0_means <- matrix(NA, 18, 6)
colnames(pmx_means) <- colnames(pmx_medians) <- colnames(pm0_means) <- 
  c("theta_y", "theta_w1","theta_w2","tau2","min_w","max_w")

# There are 18 Florida storms that we will compare
for (ste in 1:18) {
  for (pmx in 1) {
    load(paste0("rda/FL_fits/storm",ste,"_niters",niters,"krig",if(pmx==1){"pmx"},".rda"))
    
    par(mfrow=c(2,2))
    plot(fit$theta_w[,1], type="l", main=paste0("theta_w[,1] ",ste,", pmx=",pmx))
    plot(fit$theta_w[,2], type="l", main=paste0("theta_w[,2] ",ste,", pmx=",pmx))
    plot(fit$theta_y, type="l", main=paste0("theta_y ",ste,", pmx=",pmx))
    plot(fit$tau2, type="l", main=paste0("tau2 ",ste,", pmx=",pmx))
    # plot(fit$g, type="l", main=paste0("g ",ste,", pmx=",pmx))  
    
    # Save mean and median info for each TC
    if(pmx){
      pmx_means[ste,1] <- mean(fit$theta_y)
      pmx_means[ste,2:3] <- colMeans(fit$theta_w)
      pmx_means[ste,4] <- mean(fit$tau2)
      pmx_means[ste,5:6] <- range(unlist(fit$w))
      
      pmx_medians[ste,1] <- median(fit$theta_y)
      pmx_medians[ste,2:3] <- apply(fit$theta_w,2,median)
      pmx_medians[ste,4] <- median(fit$tau2)
      pmx_medians[ste,5:6] <- range(unlist(fit$w))
    } else {
      pm0_means[ste,1] <- mean(fit$theta_y)
      pm0_means[ste,2:3] <- colMeans(fit$theta_w)
      pm0_means[ste,4] <- mean(fit$tau2)
      pm0_means[ste,5:6] <- range(unlist(fit$w))
    }
  }
}

# Histograms instead of trace plots
for (ste in 1:18) {
  for (pmx in 1) {
    load(paste0("rda/FL_fits/storm",ste,"_niters",niters,"krig",if(pmx==1){"pmx"},".rda"))
    
    par(mfrow=c(2,2))
    hist(fit$theta_w[,1], main=paste0("theta_w[,1] ",ste,", pmx=",pmx))
    hist(fit$theta_w[,2], main=paste0("theta_w[,2] ",ste,", pmx=",pmx))
    hist(fit$theta_y, main=paste0("theta_y ",ste,", pmx=",pmx))
    hist(fit$tau2, main=paste0("tau2 ",ste,", pmx=",pmx))
    # hist(fit$g, main=paste0("g ",ste,", pmx=",pmx))  
  }
}

# Look at means and medians of the params across TCs
par(mfrow=c(2,3))
for (i in 1:ncol(pmx_means)) hist(pmx_means[,i], 
                                  main=paste(colnames(pmx_means)[i],"means"))
for (i in 1:ncol(pmx_medians)) hist(pmx_medians[,i], 
                                    main=paste(colnames(pmx_medians)[i], "medians"))
dev.off()
