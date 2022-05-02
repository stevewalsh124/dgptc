# modify package functions to allow true_g to be a vector

library(deepgp)
library(Rcpp)
library(parallel)
library(doParallel)

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


###################################
# Trimming function modifications #
###################################

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

#################################
# Prediction code modifications #
#################################

clean_prediction <- deepgp:::clean_prediction

predict.dgp2_SW <- function (object, x_new, lite = TRUE, store_latent = FALSE, mean_map = TRUE, 
                             EI = FALSE, cores = detectCores() - 1, precs_pred = NULL, ...){
  tic <- proc.time()[3]
  object <- clean_prediction(object)
  if(object$cov == "exp2") object$v <- 999
  if (is.numeric(x_new)) 
    x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  }
  else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) 
    warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  clusterExport(cl, c("krig_SW", "MaternFun", "eps", "invdet", "sq_dist", "Exp2Fun"))
  registerDoParallel(cl)
  thread <- NULL
  result <- foreach(thread = 1:cores) %dopar% {
    out <- list()
    if (store_latent) 
      out$w_new <- list()
    out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
    if (lite) {
      out$s2_sum <- rep(0, times = n_new)
    }
    else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) 
      out$ei_sum <- rep(0, times = n_new)
    j <- 1
    for (t in chunks[[thread]]) {
      w_t <- object$w[[t]]
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig_SW(w_t[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                       g = eps, v = object$v, precs = object$precs, precs_pred = precs_pred)
          w_new[, i] <- k$mean
        }
        else {
          k <- krig_SW(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                       g = eps, sigma = TRUE, v = object$v, precs = object$precs, precs_pred = precs_pred)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        }
      }
      if (store_latent) 
        out$w_new[[j]] <- w_new
      k <- krig_SW(object$y, sq_dist(w_t), sq_dist(w_new), 
                sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                object$tau2[t], s2 = lite, sigma = !lite, f_min = EI, 
                v = object$v, precs = object$precs, precs_pred = precs_pred)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      }
      else out$sigma_sum <- out$sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t])
        }
        else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
      j <- j + 1
    }
    return(out)
  }
  stopCluster(cl)
  mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
  if (lite) {
    s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
  }
  else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (store_latent) 
    w_new <- unlist(lapply(result, with, eval(parse(text = "w_new"))), 
                    recursive = FALSE)
  if (EI) 
    ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) 
    object$w_new <- w_new
  if (lite) {
    object$s2 <- s2_sum/object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)/precs_pred
  }
  else {
    object$Sigma <- sigma_sum/object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2)/precs_pred, n_new)
  }
  if (EI) 
    object$EI <- ei_sum/object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

krig_SW <- function (y, dx, d_new = NULL, d_cross = NULL, theta, g, tau2 = 1, 
                     s2 = FALSE, sigma = FALSE, f_min = FALSE, v = 2.5, precs, precs_pred = rep(Inf, length(y))){
  out <- list()
  if (v == 999) {
    C <- Exp2Fun(dx, c(1, theta, 0)) + diag(g/precs)
    C_cross <- Exp2Fun(d_cross, c(1, theta, 0))
  }
  else {
    C <- MaternFun(dx, c(1, theta, 0, v)) + diag(g/precs)
    C_cross <- MaternFun(d_cross, c(1, theta, 0, v))
  }
  C_inv <- invdet(C)$Mi
  out$mean <- C_cross %*% C_inv %*% y
  if (f_min) {
    if (v == 999) {
      C_cross_observed_only <- Exp2Fun(dx, c(1, theta, 0))
    }
    else C_cross_observed_only <- MaternFun(dx, c(1, theta, 0, v))
    out$f_min <- min(C_cross_observed_only %*% C_inv %*% y)
  }
  if (s2) {
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    C_new <- rep(1 + g, times = nrow(d_new))
    out$s2 <- tau2 * (C_new - diag(quadterm))
  }
  if (sigma) {
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    if (v == 999) {
      C_new <- Exp2Fun(d_new, c(1, theta, eps)) + diag(g/precs_pred) # rm diag for E(Y(X))
    }
    else C_new <- MaternFun(d_new, c(1, theta, 0, v)) + diag(g/precs_pred) # rm diag for E(Y|X)
    out$sigma <- tau2 * (C_new - quadterm)
  }
  return(out)
}


predict.dgp2vec_SW <- function (object, x_new, m = object$m, lite = TRUE, store_latent = FALSE, 
                                mean_map = TRUE, cores = detectCores() - 1, precs_pred = NULL, ...){
  tic <- proc.time()[3]
  object <- clean_prediction(object)
  if(object$cov == "exp2") object$v <- 999
  if (is.numeric(x_new)) 
    x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  if (!mean_map) 
    stop("mean_map = FALSE has not yet been implemented")
  if (lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <- NULL
  }
  else {
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m)
  }
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  }
  else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) 
    warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  clusterExport(cl, c("krig_vec_SW", "MaternFun", "eps", "invdet", "sq_dist", "create_U_SW", "MaternFun_SW",
                      "Exp2Fun", "Exp2Fun_SW"))
  registerDoParallel(cl)
  thread <- NULL
  result <- foreach(thread = 1:cores) %do% {
    out <- list()
    if (store_latent) 
      out$w_new <- list()
    out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
    if (lite) {
      out$s2_sum <- rep(0, times = n_new)
    }
    else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    j <- 1
    for (t in chunks[[thread]]) {
      w_t <- object$w[[t]]
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        k <- krig_vec_SW(w_t[, i], object$theta_w[t, i], 
                      g = eps, tau2 = 1, v = object$v, precs = object$precs, m = m, x = object$x, 
                      x_new = x_new, NNarray_pred = NN_x_new, precs_pred = precs_pred)
        w_new[, i] <- k$mean
      }
      if (store_latent) 
        out$w_new[[j]] <- w_new
      if (lite) {
        w_approx <- NULL
      }
      else {
        w_approx <- update_obs_in_approx(object$w_approx, w_t)
        w_approx <- add_pred_to_approx(w_approx, w_new, m)
      }
      k <- krig_vec_SW(object$y, object$theta_y[t], object$g[t], 
                    object$tau2[t], s2 = lite, sigma = !lite, v = object$v, precs = object$precs,
                    m = m, x = w_t, x_new = w_new, approx = w_approx, precs_pred = precs_pred)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      }
      else out$sigma_sum <- out$sigma_sum + k$sigma
      j <- j + 1
    }
    return(out)
  }
  stopCluster(cl)
  mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
  if (lite) {
    s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
  }
  else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (store_latent) 
    w_new <- unlist(lapply(result, with, eval(parse(text = "w_new"))), 
                    recursive = FALSE)
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) 
    object$w_new <- w_new
  if (lite) {
    object$s2 <- s2_sum/object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)/precs_pred
  }
  else {
    object$Sigma <- sigma_sum/object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2)/precs_pred, n_new)
  }
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

krig_vec_SW <- function (y, theta, g, tau2 = 1, s2 = FALSE, sigma = FALSE, v, precs, precs_pred = rep(Inf, length(y)), 
                         m = NULL, x = NULL, x_new = NULL, NNarray_pred = NULL, approx = NULL){
  out <- list()
  if (!sigma) {
    n_new <- nrow(x_new)
    if (is.null(NNarray_pred)) 
      NNarray_pred <- FNN::get.knnx(x, x_new, m)$nn.index
    out$mean <- vector(length = n_new)
    if (s2) 
      out$s2 <- vector(length = n_new)
    for (i in 1:n_new) {
      NN <- NNarray_pred[i, ]
      x_combined <- rbind(x[NN, , drop = FALSE], x_new[i, , drop = FALSE])
      precsub <- precs[NN]
      K <- MaternFun(sq_dist(x_combined), c(1, theta, 0, v)) + diag(c(g/precsub, g/precs_pred[i])) # rm diag for E(Y|X)
      L <- t(chol(K))
      out$mean[i] <- L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], y[NN])
      if (s2) 
        out$s2[i] <- tau2 * (L[m + 1, m + 1]^2)
    }
  }
  else {
    if (is.null(approx$observed)) 
      approx <- add_pred_to_approx(approx, x_new, m)
    yo <- y[approx$ord[approx$observed]]
    U_mat <- create_U(approx, g, theta, v, precs)
    Upp <- U_mat[!approx$observed, !approx$observed]
    Uppinv <- Matrix::solve(Upp, sparse = TRUE)
    Winv <- Matrix::crossprod(Uppinv)
    Uop <- U_mat[approx$observed, !approx$observed]
    UopTy <- Matrix::crossprod(Uop, yo)
    mu_ordered <- -Matrix::crossprod(Uppinv, UopTy)
    out$mean <- mu_ordered[approx$rev_ord_pred]
    out$sigma <- as.matrix(tau2 * Winv[approx$rev_ord_pred, 
                                       approx$rev_ord_pred])
  }
  return(out)
}
