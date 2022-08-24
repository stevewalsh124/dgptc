#########################################
# Vecchia approx function modifications #
# adapted for sample cov mtx from       #
# hierarchical model                    #
#########################################

gibbs_two_layer_vec_SW <- function (x, y, nmcmc, D, verb, initial, true_g, settings, v, 
                                    m, x_approx = NULL, w_approx = NULL, Sigma_hat) {
  print("using SW")
  if (is.null(x_approx)) 
    x_approx <- create_approx(x, m)
  if (is.null(w_approx)) 
    w_approx <- create_approx(initial$w, m)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) 
    g[1] <- initial$g
  else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w <- list()
  w[[1]] <- initial$w
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_outer <- NULL
  for (j in 2:nmcmc) {
    if (verb) 
      if (j%%1000 == 0) 
        cat(j, "\n")
    if (is.null(true_g)) {
      samp <- sample_g_vec_SW(y, g[j - 1], theta_y[j - 1], 
                              alpha = settings$alpha$g, beta = settings$beta$g, 
                              l = settings$l, u = settings$u, ll_prev = ll_outer, 
                              approx = w_approx, v = v, Sigma_hat = Sigma_hat)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else { g[j] <- true_g }
    # NOTE - outer layer unaffected by pmx
    
    samp <- sample_theta_vec_SW(y, g[j], theta_y[j - 1], alpha = settings$alpha$theta_y, 
                                beta = settings$beta$theta_y, l = settings$l, u = settings$u, 
                                outer = TRUE, ll_prev = ll_outer, approx = w_approx, 
                                v = v, tau2 = TRUE, Sigma_hat = Sigma_hat)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    # NOTE - outer layer unaffected by pmx
    
    for (i in 1:D) {
      ### NEW - sample_theta_vec_new function accepts prior mean as argument
      ### NEW - calls logl_vec_new which accepts prior mean as argument
      samp <- sample_theta_vec_new(w[[j - 1]][, i], g = eps, 
                                   theta_w[j - 1, i], alpha = settings$alpha$theta_w, 
                                   beta = settings$beta$theta_w, l = settings$l, 
                                   u = settings$u, outer = FALSE, approx = x_approx, 
                                   v = v, prior_mean = settings$w_prior_mean[, i])
      theta_w[j, i] <- samp$theta
    }
    
    samp <- sample_w_vec_SW(y, w_approx, x_approx, g[j], theta_y[j], 
                            theta_w[j, ], ll_prev = ll_outer, v = v, prior_mean = settings$w_prior_mean, Sigma_hat = Sigma_hat)
    w_approx <- samp$w_approx
    w[[j]] <- w_approx$x_ord[w_approx$rev_ord_obs, , drop = FALSE]
    ll_outer <- samp$ll
  }
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, 
              w = w, tau2 = tau2, w_approx = w_approx, x_approx = x_approx, Sigma_hat = Sigma_hat))
}

sample_g_vec_SW <- function (y, g_t, theta, alpha, beta, l, u, ll_prev = NULL, approx, v, Sigma_hat){
  g_star <- runif(1, min = l * g_t/u, max = u * g_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, approx, g_t, theta, outer = TRUE, 
                           v, Sigma_hat = Sigma_hat)$logl
  lpost_threshold <- ll_prev + dgamma(g_t - eps, alpha, beta, 
                                      log = TRUE) + log(ru) - log(g_t) + log(g_star)
  ll_new <- logl_vec_SW(y, approx, g_star, theta, outer = TRUE, 
                        v, Sigma_hat = Sigma_hat)$logl
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  }
  else {
    return(list(g = g_t, ll = ll_prev))
  }
}

sample_theta_vec_SW <- function (y, g, theta_t, alpha, beta, l, u, outer, ll_prev = NULL, 
                                 approx, v, tau2 = FALSE, Sigma_hat) {
  
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, approx, g, theta_t, outer, v, Sigma_hat = Sigma_hat)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, 
                                      log = TRUE) + log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_vec_SW(y, approx, g, theta_star, outer, v, tau2 = tau2, Sigma_hat = Sigma_hat)
  if (ll_new$logl + dgamma(theta_star, alpha, beta, log = TRUE) > 
      lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  }
  else {
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

### NEW FUNCTION
sample_theta_vec_new <- function (y, g, theta_t, alpha, beta, l, u, outer, ll_prev = NULL, 
                                  approx, v, tau2 = FALSE, prior_mean = 0) 
{
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_new(y, approx, g, theta_t, outer, v, mean = prior_mean)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t - eps, alpha, 
                                      beta, log = TRUE) + log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_vec_new(y, approx, g, theta_star, outer, v, tau2 = tau2, mean = prior_mean)
  new <- ll_new$logl + dgamma(theta_star - eps, alpha, beta, 
                              log = TRUE)
  if (new > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  }
  else {
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

### NEW FUNCTION
logl_vec_new <- function (out_vec, approx, g, theta, outer = TRUE, v, tau2 = FALSE, 
                          mean = rep(0, times = length(out_vec))) 
{
  n <- length(out_vec)
  out_vec_ord <- out_vec[approx$ord] - mean[approx$ord] ### NEW
  U_mat <- create_U(approx, g, theta, v)
  Uty <- Matrix::crossprod(U_mat, out_vec_ord)
  ytUUty <- sum(Uty^2)
  logdet <- sum(log(Matrix::diag(U_mat)))
  if (outer) {
    logl <- logdet - (n * 0.5) * log(ytUUty)
  }
  else {
    logl <- logdet - 0.5 * ytUUty
  }
  if (tau2) {
    tau2 <- c(ytUUty)/n
  }
  else tau2 <- NULL
  return(list(logl = logl, tau2 = tau2))
}

sample_w_vec_SW <- function (y, w_approx, x_approx, g, theta_y, theta_w, ll_prev = NULL, 
                             v, prior_mean, Sigma_hat) {
  
  D <- ncol(w_approx$x_ord)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, w_approx, g, theta_y, outer = TRUE, v, Sigma_hat = Sigma_hat)$logl
  for (i in 1:D) {
    w_prior <- rand_mvn_vec(x_approx, theta_w[i], v, prior_mean[, i])
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    accept <- FALSE
    count <- 0
    w_prev <- w_approx$x_ord[w_approx$rev_ord_obs, i]
    while (accept == FALSE) {
      count <- count + 1
      w_proposal <- w_prev * cos(a) + w_prior * sin(a)
      w_approx <- update_obs_in_approx(w_approx, w_proposal, 
                                       i)
      new_logl <- logl_vec_SW(y, w_approx, g, theta_y, outer = TRUE, 
                              v, Sigma_hat = Sigma_hat)$logl
      if (new_logl > ll_threshold) {
        ll_prev <- new_logl
        accept <- TRUE
      }
      else {
        if (a < 0) {
          amin <- a
        }
        else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) 
          stop("reached maximum iterations of ESS")
      }
    }
  }
  return(list(w_approx = w_approx, ll = ll_prev))
}

logl_vec_SW <- function (out_vec, approx, g, theta, outer = FALSE, v, tau2 = FALSE, Sigma_hat) {
  n <- length(out_vec)
  out_vec_ord <- out_vec[approx$ord]
  U_mat <- create_U_SW(approx, g, theta, v, Sigma_hat)
  Uty <- Matrix::crossprod(U_mat, out_vec_ord)
  ytUUty <- sum(Uty^2)
  logdet <- sum(log(Matrix::diag(U_mat)))
  if (outer) {
    logl <- logdet - (n * 0.5) * log(ytUUty)
  }
  else {
    logl <- logdet - 0.5 * ytUUty
  }
  if (tau2) {
    tau2 <- 1#c(ytUUty)/n
  }
  else tau2 <- NULL
  return(list(logl = logl, tau2 = tau2))
}

create_U_SW <- function (approx, g, theta, v, Sigma_hat) {
  n <- nrow(approx$x_ord)
  U <- U_entries_SW(approx$n_cores, n, approx$x_ord, approx$revNN,
                    approx$revCond, c(1, theta, 0, v), g = g*Sigma_hat)
  U <- c(t(U))[approx$notNA]
  U_mat <- Matrix::sparseMatrix(i = approx$col_indices, j = approx$row_pointers,
                                x = U, dims = c(n, n))
  return(U_mat)
}
