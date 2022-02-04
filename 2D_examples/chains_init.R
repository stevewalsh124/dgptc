############################################
# Create initial states for w, theta, tau2 #
# Steve Walsh                              #
# February 4 2022                          #
############################################

library(deepgp)

niters <- 25000

# save theta and tau2 (theta_y, theta_w, tau2)
init_pm0 <- init_pmx <- matrix(NA, nrow = 18, ncol = 4)
w_pm0 <- w_pmx <- list()

# There are 18 Florida storms that we will compare
for (pmx in 1) {
  for (ste in 1:18) {
    load(paste0("rda/FL_fits/storm",ste,"_niters",niters,"krig",if(pmx==1){"pmx"},".rda"))
    fit <- trim(fit, niters-1, 1) # retain only last sample
    if(pmx){
      init_pmx[ste,] <- c(fit$theta_y, fit$theta_w, fit$tau2[niters])
      w_pmx[[ste]] <- fit$w
    } else {
      init_pm0[ste,] <- c(fit$theta_y, fit$theta_w, fit$tau2[niters])
      w_pm0[[ste]] <- fit$w
    }
  }
}

init_pmx
