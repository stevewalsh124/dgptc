################################
# Compare chains b/w FL storms #
# Steve Walsh                  #
# Febraury 4 2022              #
################################

pdf("pdf/all_FL_chains.pdf")

niters <- 25000

# There are 18 Florida storms that we will compare
for (ste in 1:18) {
  for (pmx in 0:1) {
    load(paste0("rda/FL_fits/storm",ste,"_niters",niters,"krig",if(pmx==1){"pmx"},".rda"))
    
    par(mfrow=c(2,2))
    plot(fit$theta_w[,1], type="l", main=paste0("theta_w[,1] ",ste,", pmx=",pmx))
    plot(fit$theta_w[,2], type="l", main=paste0("theta_w[,2] ",ste,", pmx=",pmx))
    plot(fit$theta_y, type="l", main=paste0("theta_y ",ste,", pmx=",pmx))
    plot(fit$tau2, type="l", main=paste0("tau2 ",ste,", pmx=",pmx))
    # plot(fit$g, type="l", main=paste0("g ",ste,", pmx=",pmx))  
  }
}

dev.off()
