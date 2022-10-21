# get a matern covariance matrix for the errors
# estimates hyperpars for each (low res) run, 
# finds their averages, then builds Matern cov mtx

get_matern <- function(x, Y, nmcmc=nmcmc,nburn=nburn, v = 2.5){
  # fit one error with Matern 5/2
  e_gs <- e_tau2s <- e_thetas <- c()
  for (ri in 1:nrow(Y)) {
    e_fit <- fit_one_layer(x, Y[ri,], nmcmc = nmcmc, cov = "matern", v=v)
    e_fit <- trim(e_fit, burn = nburn, thin = kth)
    # par(mfrow=c(1,3)); plot(e_fit$g, type="l"); plot(e_fit$tau2, type="l"); plot(e_fit$theta, type="l", main = ri)
    e_gs[ri] <- mean(e_fit$g)
    e_tau2s[ri] <- mean(e_fit$tau2)
    e_thetas[ri] <- mean(e_fit$theta)
  }
  
  hist(e_gs); abline(v=mean(e_gs))
  hist(e_tau2s); abline(v=mean(e_tau2s))
  hist(sqrt(e_thetas)); abline(v=mean(e_thetas))
  
  # MaternFun
  # // distmat = matrix of SQUARED distances
  # // covparms = c(tau2, theta, g, v)
  Matern_hat <- MaternFun(plgp:::distance(x), covparms = c(mean(e_tau2s), mean(e_thetas), mean(e_gs), v))
  return(Matern_hat)
}
