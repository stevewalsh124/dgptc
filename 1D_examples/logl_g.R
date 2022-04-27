# modify package functions to allow true_g to be a vector

library(deepgp)
library(Rcpp)

sourceCpp("src/cov_SW.cpp")

##########################################
# Full likelihood function modifications #
##########################################

# link the functions from the package
check_inputs <- deepgp:::check_inputs
check_settings <- deepgp:::check_settings
check_initialization <- deepgp:::check_initialization
gibbs_two_layer <- deepgp:::gibbs_two_layer
gibbs_two_layer_vec <- deepgp:::gibbs_two_layer_vec
MaternFun <- deepgp:::MaternFun
Exp2Fun <- deepgp:::Exp2Fun
invdet <- deepgp:::invdet
sample_theta <- deepgp:::sample_theta
sample_theta_vec <- deepgp:::sample_theta_vec
eps <- sqrt(.Machine$double.eps)
create_approx <- deepgp:::create_approx
rand_mvn_vec <- deepgp:::rand_mvn_vec
update_obs_in_approx <- deepgp:::update_obs_in_approx


fit_two_layer_SW <- function (x, y, D = ifelse(is.matrix(x), ncol(x), 1), nmcmc = 10000, 
          verb = TRUE, w_0 = NULL, g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
          true_g = NULL, settings = NULL, cov = c("matern", "exp2"), 
          v = 2.5, vecchia = FALSE, m = min(25, length(y) - 1), precs = NULL) {
  tic <- proc.time()[3]
  cov <- match.arg(cov)
  if (vecchia & cov == "exp2") {
    message("vecchia = TRUE requires matern covariance, proceeding with cov = 'matern'")
    cov <- "matern"
  }
  if (vecchia & sum(duplicated(x)) > 0) 
    stop("unable to handle duplicate training locations")
  if (is.numeric(x)) 
    x <- as.matrix(x)
  test <- check_inputs(x, y, true_g)
  settings <- check_settings(settings, layers = 2, D, length(y))
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 2, x = x, 
                                  D = D, vecchia = vecchia, v = v, m = m)
  if (m >= length(y)) 
    stop("m must be less than the length of y")
  if (cov == "matern") 
    if (!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, 
              cov = cov)
  if (cov == "matern") 
    out$v <- v
  if (vecchia) 
    out$m <- m
  if (vecchia) {
    if(is.null(precs)){
      samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, initial, 
                                 true_g, settings, v, m)
    } else {
      if(length(precs)==nrow(x)){
        samples <- gibbs_two_layer_vec_SW(x, y, nmcmc, D, verb, initial, 
                                      true_g, settings, v, m, precs = precs)
      } else {print(stop("length of precs must be length/nrow of x"))}
    }
  } else {
    # check if length(precs) is NULL, or nrow(as.matrix(x))
    if(is.null(precs)){
      samples <- gibbs_two_layer(x, y, nmcmc, D, verb, initial, 
                                 true_g, settings, cov, v)
    } else {
      if(length(precs)==nrow(x)){
        samples <- gibbs_two_layer_SW(x, y, nmcmc, D, verb, initial, 
                                   true_g, settings, cov, v, precs)
      } else {print(stop("length of precs must be the length/nrow of x"))}
    }
  }
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  if(!is.null(precs)) 
    out$precs <- precs
  if (vecchia) 
    class(out) <- "dgp2vec"
  else class(out) <- "dgp2"
  return(out)
}


# use this version when g is a vector
gibbs_two_layer_SW <- function (x, y, nmcmc, D, verb, initial, true_g, settings, cov, v, precs){
  print("using SW")
  dx <- sq_dist(x)
  dw <- sq_dist(initial$w)
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
      samp <- sample_g_SW(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v, cov = cov, precs = precs)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else {g[j] <- true_g}
    samp <- sample_theta_SW(y, dw, g[j], theta_y[j - 1], alpha = settings$alpha$theta_y, 
                         beta = settings$beta$theta_y, l = settings$l, u = settings$u, 
                         outer = TRUE, ll_prev = ll_outer, v = v, cov = cov, 
                         tau2 = TRUE, precs = precs)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) 
      tau2[j] <- tau2[j - 1]
    else tau2[j] <- samp$tau2
    for (i in 1:D) {
      samp <- sample_theta(w[[j - 1]][, i], dx, g = eps, 
                           theta_w[j - 1, i], alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v)
      theta_w[j, i] <- samp$theta
    }
    samp <- sample_w_SW(y, w[[j - 1]], dw, dx, g[j], theta_y[j], 
                     theta_w[j, ], ll_prev = ll_outer, v = v, cov = cov, 
                     prior_mean = settings$w_prior_mean, precs = precs)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  }
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, 
              w = w, tau2 = tau2))
}

sample_g_SW <- function (out_vec, in_dmat, g_t, theta, alpha, beta, l, u, ll_prev = NULL, v, cov, precs){
  g_star <- runif(1, min = l * g_t/u, max = u * g_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, in_dmat, g_t, theta, outer = TRUE, 
                    v, cov, precs = precs)$logl
  lpost_threshold <- ll_prev + dgamma(g_t - eps, alpha, beta, 
                                      log = TRUE) + log(ru) - log(g_t) + log(g_star)
  ll_new <- logl_SW(out_vec, in_dmat, g_star, theta, outer = TRUE, 
                 v, cov, precs = precs)$logl
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  }
  else {
    return(list(g = g_t, ll = ll_prev))
  }
}

# change sample_theta for the outer=TRUE case
sample_theta_SW <- function (out_vec, in_dmat, g, theta_t, alpha, beta, l, u, outer, 
                                   ll_prev = NULL, v, cov, tau2 = FALSE, precs){
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, in_dmat, g, theta_t, outer, v, cov, precs = precs)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, 
                                      log = TRUE) + log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_SW(out_vec, in_dmat, g, theta_star, outer, v, 
                 cov, tau2 = tau2, precs = precs)
  if (ll_new$logl + dgamma(theta_star, alpha, beta, log = TRUE) > 
      lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  }
  else {
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

# change sample_w (not the prior, just logl)
sample_w_SW <- function (out_vec, w_t, w_t_dmat, in_dmat, g, theta_y, theta_w, 
          ll_prev = NULL, v, cov, prior_mean, precs) {
  D <- ncol(w_t)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, w_t_dmat, g, theta_y, outer = TRUE, 
                    v = v, cov = cov, precs = precs)$logl
  count <- vector(length = D)
  for (i in 1:D) {
    if (cov == "matern") {
      w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], sigma = MaternFun(in_dmat, c(1, theta_w[i], 0, v)))
    }
    else w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], sigma = Exp2Fun(in_dmat, c(1, theta_w[i], 0)))
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    accept <- FALSE
    count <- 0
    w_prev <- w_t[, i]
    while (accept == FALSE) {
      count <- count + 1
      w_t[, i] <- w_prev * cos(a) + w_prior * sin(a)
      dw <- sq_dist(w_t)
      new_logl <- logl_SW(out_vec, dw, g, theta_y, outer = TRUE, 
                       v = v, cov = cov, precs = precs)$logl
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
  return(list(w = w_t, ll = ll_prev, dw = dw))
}

# change the log likelihood evaluation
# MaternFun, Exp2Fun: didn't change these, just add the diagonal outside of C code
logl_SW <- function (out_vec, in_dmat, g, theta, outer = TRUE, v, cov, tau2 = FALSE, precs){
    n <- length(out_vec)
    if (cov == "matern") {
      K <- MaternFun(in_dmat, c(1, theta, 0, v)) + diag(g/precs) #tau2=1 here, nugget is tau2*g
    }
    else K <- Exp2Fun(in_dmat, c(1, theta, eps)) + diag(g/precs) #tau2=1 here, nugget is tau2*g
    id <- invdet(K)
    quadterm <- t(out_vec) %*% id$Mi %*% (out_vec)
    if (outer) {
      logl <- (-n * 0.5) * log(quadterm) - 0.5 * id$ldet
    }
    else logl <- (-0.5) * id$ldet - 0.5 * quadterm
    if (tau2) {
      tau2 <- c(quadterm)/n
    }
    else tau2 <- NULL
    return(list(logl = c(logl), tau2 = tau2))
}


#########################################
# Vecchia approx function modifications #
#########################################

gibbs_two_layer_vec_SW <- function (x, y, nmcmc, D, verb, initial, true_g, settings, v, 
          m, x_approx = NULL, w_approx = NULL, precs) {
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
                           approx = w_approx, v = v, precs = precs)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else { g[j] <- true_g }
    samp <- sample_theta_vec_SW(y, g[j], theta_y[j - 1], alpha = settings$alpha$theta_y, 
                             beta = settings$beta$theta_y, l = settings$l, u = settings$u, 
                             outer = TRUE, ll_prev = ll_outer, approx = w_approx, 
                             v = v, tau2 = TRUE, precs = precs)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    for (i in 1:D) {
      samp <- sample_theta_vec(w[[j - 1]][, i], g = eps, 
                               theta_w[j - 1, i], alpha = settings$alpha$theta_w, 
                               beta = settings$beta$theta_w, l = settings$l, 
                               u = settings$u, outer = FALSE, approx = x_approx, 
                               v = v)
      theta_w[j, i] <- samp$theta
    }
    samp <- sample_w_vec_SW(y, w_approx, x_approx, g[j], theta_y[j], 
                         theta_w[j, ], ll_prev = ll_outer, v = v, prior_mean = settings$w_prior_mean, precs = precs)
    w_approx <- samp$w_approx
    w[[j]] <- w_approx$x_ord[w_approx$rev_ord_obs, , drop = FALSE]
    ll_outer <- samp$ll
  }
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, 
              w = w, tau2 = tau2, w_approx = w_approx, x_approx = x_approx, precs = precs))
}

sample_g_vec_SW <- function (y, g_t, theta, alpha, beta, l, u, ll_prev = NULL, approx, v, precs){
  g_star <- runif(1, min = l * g_t/u, max = u * g_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, approx, g_t, theta, outer = TRUE, 
                        v, precs = precs)$logl
  lpost_threshold <- ll_prev + dgamma(g_t - eps, alpha, beta, 
                                      log = TRUE) + log(ru) - log(g_t) + log(g_star)
  ll_new <- logl_vec_SW(y, approx, g_star, theta, outer = TRUE, 
                     v, precs = precs)$logl
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  }
  else {
    return(list(g = g_t, ll = ll_prev))
  }
}

sample_theta_vec_SW <- function (y, g, theta_t, alpha, beta, l, u, outer, ll_prev = NULL, 
          approx, v, tau2 = FALSE, precs) {

  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, approx, g, theta_t, outer, v, precs = precs)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, 
                                      log = TRUE) + log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_vec_SW(y, approx, g, theta_star, outer, v, tau2 = tau2, precs = precs)
  if (ll_new$logl + dgamma(theta_star, alpha, beta, log = TRUE) > 
      lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  }
  else {
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

sample_w_vec_SW <- function (y, w_approx, x_approx, g, theta_y, theta_w, ll_prev = NULL, 
          v, prior_mean, precs) {

  D <- ncol(w_approx$x_ord)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec_SW(y, w_approx, g, theta_y, outer = TRUE, v, precs = precs)$logl
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
                           v, precs = precs)$logl
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

logl_vec_SW <- function (out_vec, approx, g, theta, outer = TRUE, v, tau2 = FALSE, precs) {
  n <- length(out_vec)
  out_vec_ord <- out_vec[approx$ord]
  U_mat <- create_U_SW(approx, g, theta, v, precs)
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

create_U_SW <- function (approx, g, theta, v, precs) {
  n <- nrow(approx$x_ord)
  U <- U_entries_SW(approx$n_cores, n, approx$x_ord, approx$revNN,
                 approx$revCond, c(1, theta, 0, v), g = g/precs)
  U <- c(t(U))[approx$notNA]
  U_mat <- Matrix::sparseMatrix(i = approx$col_indices, j = approx$row_pointers,
                                x = U, dims = c(n, n))
  return(U_mat)
}


###########################################
# Post-processing functions modifications #
###########################################

# modify trim.dgp2; trim.dgp2vec <- trim.dgp2
trim_SW <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  if(is.matrix(object$g)){  
    object$g <- object$g[indx,, drop = FALSE]
  } else {
    object$g <- object$g[indx, drop = FALSE]
  }
  object$theta_y <- object$theta_y[indx, drop = FALSE]
  object$theta_w <- object$theta_w[indx, , drop = FALSE]
  object$w <- as.list(object$w[indx])
  object$tau2 <- object$tau2[indx, drop = FALSE]
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

