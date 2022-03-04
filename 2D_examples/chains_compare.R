################################
# Compare chains b/w FL storms #
# Steve Walsh                  #
# Febraury 4 2022              #
################################

library(deepgp)

# Change load() call to remove /orig_50k if necessary
# Change pdf name for each batch

# Choose based on the files you want to load
# Save every kth iteration (thinning MCMC)
niters <- 25000
k <- 100

pdf(paste0("pdf/all_FL_chains_",niters,"_thin",k,".pdf"))

pmx_means <- pmx_medians <- pm0_means <- matrix(NA, 18, 6)
colnames(pmx_means) <- colnames(pmx_medians) <- colnames(pm0_means) <- 
  c("theta_y", "theta_w1","theta_w2","tau2","min_w","max_w","g")

# There are 18 Florida storms that we will compare
# Trace plots, histograms and acf plots of the chains 
# for theta = (theta_y, theta_w1, theta_w2) and tau2
for (ste in 1:18) {
  print(ste)
  for (pmx in 1) {
    load(paste0("rda/FL_fits/storm",ste,"_niters",niters,".rda"))
    fit <- trim(fit, 0, k)

    par(mfrow=c(2,3))
    plot(fit$theta_y, type="l", main=paste0("theta_y ",ste,", pmx=",pmx))
    plot(fit$theta_w[,1], type="l", main=paste0("theta_w[,1] ",ste,", pmx=",pmx))
    plot(fit$theta_w[,2], type="l", main=paste0("theta_w[,2] ",ste,", pmx=",pmx))
    plot(fit$tau2, type="l", main=paste0("tau2 ",ste,", pmx=",pmx))
    plot(fit$g, type="l", main=paste0("g ",ste,", pmx=",pmx))
    
    par(mfrow=c(2,3))
    hist(fit$theta_y, main=paste0("theta_y ",ste,", pmx=",pmx))
    hist(fit$theta_w[,1], main=paste0("theta_w[,1] ",ste,", pmx=",pmx))
    hist(fit$theta_w[,2], main=paste0("theta_w[,2] ",ste,", pmx=",pmx))
    hist(fit$tau2, main=paste0("tau2 ",ste,", pmx=",pmx))
    hist(fit$g, main=paste0("g ",ste,", pmx=",pmx))
    
    par(mfrow=c(2,3))
    acf(fit$theta_y, lag.max = length(fit$theta_y))
    acf(fit$theta_w[,1], lag.max = length(fit$theta_w[,1]))
    acf(fit$theta_w[,2], lag.max = length(fit$theta_w[,2]))
    acf(fit$tau2, lag.max = length(fit$tau2))
    acf(fit$g, lag.max = length(fit$g))
    
    # Save mean and median info for each TC,
    # as well as min and max over all w's per TC
    if(pmx){
      pmx_means[ste,1] <- mean(fit$theta_y)
      pmx_means[ste,2:3] <- colMeans(fit$theta_w)
      pmx_means[ste,4] <- mean(fit$tau2)
      pmx_means[ste,5:6] <- range(unlist(fit$w))
      pmx_means[ste,7] <- mean(fit$g)
      
      pmx_medians[ste,1] <- median(fit$theta_y)
      pmx_medians[ste,2:3] <- apply(fit$theta_w,2,median)
      pmx_medians[ste,4] <- median(fit$tau2)
      pmx_medians[ste,5:6] <- range(unlist(fit$w))
      pmx_means[ste,7] <- median(fit$g)
      
    } else {
      pm0_means[ste,1] <- mean(fit$theta_y)
      pm0_means[ste,2:3] <- colMeans(fit$theta_w)
      pm0_means[ste,4] <- mean(fit$tau2)
      pm0_means[ste,5:6] <- range(unlist(fit$w))
      pm0_means[ste,7] <- mean(fit$g)
      
    }
  }
}

# Save 
save(pmx_means, file = paste0("rda/FL_summaries/pmx_means_",niters,"_thin",k,".rda"))
save(pmx_medians, file = paste0("rda/FL_summaries/pmx_medians_",niters,"_thin",k,".rda"))

# Look at means and medians of the params across TCs
par(mfrow=c(2,3))
for (i in (1:ncol(pmx_means))[-5]) hist(pmx_means[,i], 
                                  main=paste(colnames(pmx_means)[i],"means"))
for (i in (1:ncol(pmx_medians)[-6])) hist(pmx_medians[,i], 
                                    main=paste(colnames(pmx_medians)[i], "medians"))
dev.off()
