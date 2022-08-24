###################################
# Prediction code modifications   #
# adapted for sample cov mtx from #
# hierarchical model              #
###################################

clean_prediction <- deepgp:::clean_prediction

predict.dgp2_SW <- function (object, x_new, lite = TRUE, store_latent = FALSE, mean_map = TRUE, 
                             EI = FALSE, cores = detectCores() - 1, precs_pred = NULL, ...){
  tic <- proc.time()[3]
  # print(object$cov); print(object$v); print(object$nmcmc)
  # object <- clean_prediction(object)
  # print(object$cov); print(object$v); print(object$nmcmc)
  if(object$cov == "exp2"){object$v <- 999}
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
                       g = eps, v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
          w_new[, i] <- k$mean
        }
        else {
          k <- krig_SW(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                       g = eps, sigma = TRUE, v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        }
      }
      if (store_latent) 
        out$w_new[[j]] <- w_new
      k <- krig_SW(object$y, sq_dist(w_t), sq_dist(w_new), 
                   sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                   object$tau2[t], s2 = lite, sigma = !lite, f_min = EI, 
                   v = object$v, Sigma_hat = object$Sigma_hat, precs_pred = precs_pred)
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
                     s2 = FALSE, sigma = FALSE, f_min = FALSE, v = 2.5, Sigma_hat, precs_pred = rep(Inf, length(y))){
  out <- list()
  n <- length(y)
  if (v == 999) {
    C <- Exp2Fun(dx, c(1, theta, g)) + Sigma_hat/tau2 #+ diag(x = eps, nrow = n)
    C_cross <- Exp2Fun(d_cross, c(1, theta, 0))
  }
  else {
    C <- MaternFun(dx, c(1, theta, 0, v)) + Sigma_hat/tau2 #+ diag(x = eps, nrow = n)
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
      C_new <- Exp2Fun(d_new, c(1, theta, g)) + diag(1/precs_pred)/tau2 # rm diag for E(Y(X))
    }
    else C_new <- MaternFun(d_new, c(1, theta, g, v)) + diag(1/precs_pred)/tau2 # rm diag for E(Y|X)
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
  result <- foreach(thread = 1:cores) %dopar% {
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
