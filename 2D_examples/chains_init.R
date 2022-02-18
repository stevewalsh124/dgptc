############################################
# Create initial states for w, theta, tau2 #
# Steve Walsh                              #
# February 4 2022                          #
############################################

library(deepgp)

# Choose based on the files you want to load
niters <- 25000

# Gather param estimates when prior mean for W is x (T) or 0 (F)
pmx <- T

# save theta and tau2 (theta_y, theta_w, tau2)
init_param <- matrix(NA, nrow = 18, ncol = 4)
init_w <- list()

# There are 18 Florida storms that we will compare
for (ste in 1:18) {
  load(paste0("rda/FL_fits/storm",ste,"_niters",niters,".rda"))
  fit <- trim(fit, niters-1, 1) # retain only last sample
  init_param[ste,] <- c(fit$theta_y, fit$theta_w, fit$tau2)
  init_w[[ste]] <- fit$w
}

colnames(init_param) <- c("theta_y","theta_w1","theta_w2","tau2")

# Save each storm's last iteration of param values
save(init_param, file = paste0("rda/burn_params_FL",if(pmx==1){"pmx"},".rda"))
save(init_w, file = paste0("rda/burn_w_FL",if(pmx==1){"pmx"},".rda"))
