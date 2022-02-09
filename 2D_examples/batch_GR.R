#################################
# Compare chains (Gelman-Rubin) #
# Steve Walsh                   #
# February 9 2022               #
#################################

library(Bolstad2) #GelmanRubin
library(coda) # gelman.diag

GRs <- gds <- matrix(NA, 18, 4)

for (ste in c(2,4,6)) {
  
  # First 25k (burn-in)
  load(paste0("rda/FL_fits/storm",ste,"_niters25000krigpmx.rda"))
  fit_burn <- fit
  
  # First 50k (after burn-in)
  load(paste0("rda/FL_fits/orig_50k/storm",ste,"_niters50000krigpmx.rda"))
  fit_orig <- fit
  
  # Combine the 25 burn-in with 50k samples after
  l <- list(fit_burn, fit_orig)
  keys <- unique(unlist(lapply(l, names)))
  fit_orig <- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)
  fit_orig$theta_w <- matrix(fit_orig$theta_w, 75000, 2)
  
  # # Next 50k (100k samples after "burn-in")
  # load(paste0("rda/FL_fits/storm",ste,"_niters50000krigpmx.rda"))
  # fit_newer <- fit
  
  # Compare with 75k chain
  load(paste0("rda/FL_fits/storm",ste,"_niters75000krigpmx.rda"))
  fit_newer <- fit
  
  rm(fit)
  
  # Check to make sure that the last of the prev chain is the 1st of the next chain
  # fit_burn$theta_w[25000,]
  # fit_orig$theta_w[1,]
  # fit_orig$theta_w[50000,]
  # fit_newer$theta_w[1,]
  # fit_newer$theta_w[50000,]
  
  GRs[ste,1] <- GelmanRubin(cbind(fit_orig$theta_y, fit_newer$theta_y))[[5]]
  GRs[ste,2] <- GelmanRubin(cbind(fit_orig$theta_w[,1], fit_newer$theta_w[,1]))[[5]]
  GRs[ste,3] <- GelmanRubin(cbind(fit_orig$theta_w[,2], fit_newer$theta_w[,2]))[[5]]
  GRs[ste,4] <- GelmanRubin(cbind(fit_orig$tau2, fit_newer$tau2))[[5]]
  
  gds[ste,1] <- gelman.diag(mcmc(data=mcmc.list(mcmc(fit_orig$theta_y), 
                                                mcmc(fit_newer$theta_y))))$psrf[1]
  gds[ste,2] <- gelman.diag(mcmc(data=mcmc.list(mcmc(fit_orig$theta_w[,1]), 
                                                mcmc(fit_newer$theta_w[,1]))))$psrf[1]
  gds[ste,3] <- gelman.diag(mcmc(data=mcmc.list(mcmc(fit_orig$theta_w[,2]), 
                                                mcmc(fit_newer$theta_w[,2]))))$psrf[1]
  gds[ste,4] <- gelman.diag(mcmc(data=mcmc.list(mcmc(fit_orig$tau2), 
                                                mcmc(fit_newer$tau2))))$psrf[1]
}

GRs[complete.cases(GRs),]

gds[complete.cases(gds),]
