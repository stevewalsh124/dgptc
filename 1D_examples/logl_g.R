# modify package functions to allow true_g to be a vector

library(deepgp)

bte <- 3 # cols 3-18 are low res

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

step <- 499
i <- 1 # Model 1, choose from 000-111
pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
vars <- 1/prec_lowres[index_list$lowres.ix]

# wavenumber is X, a particular lowres run in Y
x <- log10(k[index_list$lowres.ix])
y <- log10(pk2[index_list$lowres.ix, bte])
x <- (x - min(x))/(max(x)-min(x))
y <- (y - mean(y))/sd(y)

check_inputs <- deepgp:::check_inputs
check_settings <- deepgp:::check_settings
check_initialization <- deepgp:::check_initialization
gibbs_two_layer <- deepgp:::gibbs_two_layer
MaternFun <- deepgp:::MaternFun
Exp2Fun <- deepgp:::Exp2Fun
invdet <- deepgp:::invdet
sample_theta <- deepgp:::sample_theta
eps <- sqrt(.Machine$double.eps)

fit_two_layer_SW <- function (x, y, D = ifelse(is.matrix(x), ncol(x), 1), nmcmc = 10000, 
          verb = TRUE, w_0 = NULL, g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
          true_g = NULL, settings = NULL, cov = c("matern", "exp2"), 
          v = 2.5, vecchia = FALSE, m = min(25, length(y) - 1)) {
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
                                  D = D, vecchia = vecchia, cov = cov, v = v, m = m)
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
    samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, 
                                   initial, true_g, settings, v, m)
  }
  else {
    # check if length(true_g) is equal to 1, or nrow(as.matrix(x))
    if(length(true_g)==1){
      samples <- gibbs_two_layer(x, y, nmcmc, D, verb, initial, 
                                 true_g, settings, cov, v)
    } else {
      if(length(true_g)==nrow(x)){
        samples <- gibbs_two_layer_SW(x, y, nmcmc, D, verb, initial, 
                                   true_g, settings, cov, v)
      } else {print(stop("length of true_g must be 1 or length/nrow of x"))}
    }
  }
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  if (vecchia) 
    class(out) <- "dgp2vec"
  else class(out) <- "dgp2"
  return(out)
}

# use this version when g is a vector
gibbs_two_layer_SW <- function (x, y, nmcmc, D, verb, initial, true_g, settings, cov, v){
  print("using SW")
  dx <- sq_dist(x)
  dw <- sq_dist(initial$w)
  if (is.null(true_g)) {g[1] <- initial$g} else {
        if(length(true_g) != nrow(x)) stop("true_g should be length 1 or length of x")
        g <- matrix(nrow = nmcmc, ncol = length(true_g))
        g[1,] <- true_g
  }
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
      if (j%%500 == 0) 
        cat(j, "\n")
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v, cov = cov)
      g[j] <- samp$g
      ll_outer <- samp$ll
    }
    else {
      if(length(true_g)==1) {g[j] <- true_g} else {g[j,] <- true_g}
    }
    samp <- sample_theta_SW(y, dw, g[j,], theta_y[j - 1], alpha = settings$alpha$theta_y, 
                         beta = settings$beta$theta_y, l = settings$l, u = settings$u, 
                         outer = TRUE, ll_prev = ll_outer, v = v, cov = cov, 
                         tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) 
      tau2[j] <- tau2[j - 1]
    else tau2[j] <- samp$tau2
    for (i in 1:D) {
      samp <- sample_theta(w[[j - 1]][, i], dx, g = eps, 
                           theta_w[j - 1, i], alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v, cov = cov)
      theta_w[j, i] <- samp$theta
    }
    samp <- sample_w_SW(y, w[[j - 1]], dw, dx, g[j,], theta_y[j], 
                     theta_w[j, ], ll_prev = ll_outer, v = v, cov = cov, 
                     prior_mean = settings$w_prior_mean)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  }
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, 
              w = w, tau2 = tau2))
}

# change sample_theta for the outer=TRUE case
sample_theta_SW <- function (out_vec, in_dmat, g, theta_t, alpha, beta, l, u, outer, 
                                   ll_prev = NULL, v, cov, tau2 = FALSE){
  theta_star <- runif(1, min = l * theta_t/u, max = u * theta_t/l)
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, in_dmat, g, theta_t, outer, v, cov)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, 
                                      log = TRUE) + log(ru) - log(theta_t) + log(theta_star)
  ll_new <- logl_SW(out_vec, in_dmat, g, theta_star, outer, v, 
                 cov, tau2 = tau2)
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
          ll_prev = NULL, v, cov, prior_mean) {
  D <- ncol(w_t)
  if (is.null(ll_prev)) 
    ll_prev <- logl_SW(out_vec, w_t_dmat, g, theta_y, outer = TRUE, 
                    v = v, cov = cov)$logl
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
                       v = v, cov = cov)$logl
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
logl_SW <- function (out_vec, in_dmat, g, theta, outer = TRUE, v, cov, tau2 = FALSE){
    n <- length(out_vec)
    if (cov == "matern") {
      K <- MaternFun(in_dmat, c(1, theta, 0, v)) + diag(g)
    }
    else K <- Exp2Fun(in_dmat, c(1, theta, eps)) + diag(g)
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

# MaternFun, Exp2Fun: don't change these, just add the diagonal outside of C code

# try it out

# check to make sure results are the same for scalar g
set.seed(1)
fit <- fit_two_layer(x, y, nmcmc = 1000, true_g = 1e-4)
set.seed(1)
fit2 <- fit_two_layer_SW(x, y, nmcmc = 1000, true_g = rep(1e-4, nrow(as.matrix(x))))

all.equal(fit$theta_y, fit2$theta_y)
