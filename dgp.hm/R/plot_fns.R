library(scoringRules)

plot.warp <- function(fit, wl = 1, wh = length(fit$x), ref.scale = 1){
  x <- fit$x
  
  for (j in 1:fit$nmcmc) {
    y <- c(fit$w[[j]])
    
    # Second: undo scale
    # find the length of the deformed reference vector
    y.ref <- y[wl] - y[wh]
    scale.y <- sqrt(t(y.ref)%*%y.ref)
    # numerator comes from length of vector with points zz, zo
    scale.back =  c(ref.scale/scale.y)#1/sqrt(sum(y.translate[,2]^2))
    y.rot.scale = y * scale.back
    
    # Third: Translate back
    translation <- y.rot.scale[wl] - x[wl]
    y.rot.scale.tran <- y.rot.scale - translation
    wstar <- y.rot.scale.tran
    
    # flip the graph if necessary in either direction
    # remove the reflections
    if(wstar[wl]>wstar[wh]) {wstar <- (wstar-wstar[wl])/(wstar[wh]-wstar[wl])}
    fit$w[[j]] <- wstar
  }
  
  ws <- do.call(rbind, fit$w)
  w_avg <- colMeans(ws)
  w_lb <- apply(ws, 2, function(x) quantile(x, 0.025))
  w_ub <- apply(ws, 2, function(x) quantile(x, 0.975))
  
  par(mfrow=c(1,1))
  plot(x, w_avg, type="l", col="blue", main = "warpings, with mean and 95% CI")
  # for(j in 1:fit$nmcmc) lines(x, fit$w[[j]])
  # lines(x, w_avg, col="blue")
  lines(x, w_lb, col="blue", lty=2)
  lines(x, w_ub, col="blue", lty=2)
  abline(0, 1, col="red")
}

plot.krig <- function(fit, zz=fit$mean, Y=parent.frame()$Y, precs_pred=parent.frame()$precs_pred){
  
  rg <- range(c(zz-2*sqrt(fit$s2_smooth*mean(fit$tau2)*mean(fit$g)),
                zz+2*sqrt(fit$s2_smooth*mean(fit$tau2)*mean(fit$g)),
                zz-2*sqrt(fit$s2*mean(fit$tau2)),
                zz+2*sqrt(fit$s2*mean(fit$tau2)),
                zz-2*sqrt(1/precs_pred),
                zz+2*sqrt(1/precs_pred), Y))
  
  plot(fit$x_new, zz, type="l", col="blue", lwd=1.5,
       ylim = rg, main = "fit")
  for (i in 1:nrow(Y))  lines(c(x), Y[i,], col="gray")#-y_avg)
  lines(fit$x_new, zz, col="blue", lwd=1.5)
  lines(fit$x_new, zz-2*sqrt(fit$s2_smooth), col=2, lwd=1.5)
  lines(fit$x_new, zz+2*sqrt(fit$s2_smooth), col=2, lwd=1.5)
  lines(fit$x_new, zz-2*sqrt(fit$s2), col=2, lwd=1.5)
  lines(fit$x_new, zz+2*sqrt(fit$s2), col=2, lwd=1.5)
  lines(fit$x_new, zz-2*sqrt(1/precs_pred), col=3, lwd=1.5)
  lines(fit$x_new, zz+2*sqrt(1/precs_pred), col=3, lwd=1.5)
  lines(fit$x_new, zz-2*sqrt(1/(nrun*precs_pred)), col=3, lwd=1.5)
  lines(fit$x_new, zz+2*sqrt(1/(nrun*precs_pred)), col=3, lwd=1.5)
  if(exists("y_hi")) lines(fit$x, y_hi, lwd=1.5, lty=2)
  if(exists("ytrue")) lines(fit$x, ytrue, lwd=1.5, lty=2)
  legend("topright", col=c("blue",2:3,1), legend = c("mean","UQ", "precs (low res)","true S"), 
         lty=c(1,1,1,2), lwd=1.5)
}

# Plot simulations of truth given the data
plot.true <- function(fit, S_e = fit$Sigma_hat, ne = 1, tolpower = -10,
                      Y = parent.frame()$Y, nrun = nrow(Y)){
  S_ei <- matrix.Moore.Penrose2(S_e, tolp = tolpower)
  S_ei <- (S_ei+t(S_ei))/2
  
  Cs <- matrix(NA, length(fit$x)^2, fit$nmcmc)
  Ms <- matrix(NA, length(fit$x), fit$nmcmc)
  Ss <- St <- Sw <- Sx <- matrix(NA, length(fit$x), fit$nmcmc*ne)
  for (i in 1:fit$nmcmc){
    if( i %% 2500 == 0) print(i)
    theta <- fit$theta_y[i]
    # tau2 <- fit$tau2[i]
    # t2S_ei <- 1/tau2 * S_ei
    g <- fit$g[i]
    W <- fit$w[[i]]
    dw <- sq_dist(W)
    if(v==999){
      S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
    } else {
      S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
    }
    S_si <- (S_si+t(S_si))/2
    C <- Cs[,i] <- solve(S_si + S_ei)
    C <- (C+t(C))/2
    M <- Ms[,i] <- C %*% S_ei %*% fit$y
    if(i %% 1000 ==0) print(range(M))
    # Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fit$y))
    Ss[,((i-1)*ne+1):(i*ne)] <- t(mvtnorm::rmvnorm(n = ne, mean = M, sigma = C, method = "eigen"))
    # Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
    # Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
  }
  
  m <- rowMeans(Ss) #- y_avg
  ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
  lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
  ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
  lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg
  
  if(exists("seed") & exists("ytrue")) write.csv(logs_sample(y = ytrue, dat = Ss), file = paste0("csv/logS/notap/logscore_",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  if(exists("y_hi")) emp_cover <- round(mean(y_hi > lb & y_hi < ub),3)
  if(exists("y_hi")) emp_cover99 <- round(mean(y_hi > lbb & y_hi < ubb),3)
  
  if(exists("ytrue")) emp_cover <- round(mean(ytrue > lb & ytrue < ub),3)
  if(exists("ytrue")) emp_cover99 <- round(mean(ytrue > lbb & ytrue < ubb),3)
  if(exists("seed") & exists("ytrue")) write.csv(mean((m-ytrue)^2), file = paste0("csv/MSE/notap/",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  if(exists("seed") & exists("ytrue")) write.csv(c(emp_cover, emp_cover99), file = paste0("csv/emp_cover_",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  plot(fit$x, fit$y, type="n",
       ylim = range(c(m, lb, ub, lbb, ubb, Y)), 
       main = paste0("est both,", if(exists("emp_cover")){emp_cover}, " ", if(exists("emp_cover99")){emp_cover99}))
  
  for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  if(exists("y_hi")) lines(fit$x, y_hi, lwd=1.5, col="red")
  if(exists("ytrue")) lines(fit$x, ytrue, lwd=1.5, col="red")
  lines(fit$x, fit$y, lwd=1.5, lty=2)
  lines(fit$x, m , col="blue")
  lines(fit$x, lb , col="blue", lty=2)
  lines(fit$x, ub , col="blue", lty=2)
  lines(fit$x, lbb , col="darkblue", lty=2)
  lines(fit$x, ubb , col="darkblue", lty=2)
  legend("bottomright", legend = c("true", "sample avg", "UQ"), 
         col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
  # Using law of total expectation and variance instead
  mu_bar <- rowMeans(Ms)
  C_bar <- matrix(rowMeans(Cs),  length(fit$x),  length(fit$x))
  cov_mut <- cov(t(Ms))
  
  plot(x, mu_bar, type="n", main = "law of total E & V; MPI version", col="blue",
       ylim=range(Y))
  for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  lines(x, mu_bar, col="blue")
  if(exists("y_hi")) lines(x, y_hi, col="red")
  if(exists("ytrue")) lines(x, ytrue, col="red")
  lines(x, fit$y, lty=2)
  lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  legend("bottomright", legend = c("true", "sample avg", "UQ"), 
         col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
}


# Plot simulations of truth given the data
plot.true.combo <- function(fit, S_e = fit$Sigma_hat, ne = 1, tolpower = -10,
                      Y = parent.frame()$Y, nrun = nrow(Y)){
  S_ei <- matrix.Moore.Penrose2(S_e, tolp = tolpower)
  S_ei <- (S_ei+t(S_ei))/2
  
  Cs <- matrix(NA, length(fit$x)^2, fit$nmcmc)
  Ms <- matrix(NA, length(fit$x), fit$nmcmc)
  Ss <- St <- Sw <- Sx <- matrix(NA, length(fit$x), fit$nmcmc*ne)
  for (i in 1:fit$nmcmc){
    if( i %% 2500 == 0) print(i)
    theta <- fit$theta_y[i]
    # tau2 <- fit$tau2[i]
    # t2S_ei <- 1/tau2 * S_ei
    g <- fit$g[i]
    W <- fit$w[[i]]
    dw <- sq_dist(W)
    if(v==999){
      S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
    } else {
      S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
    }
    S_si <- (S_si+t(S_si))/2
    C <- Cs[,i] <- solve(S_si + S_ei)
    C <- (C+t(C))/2
    M <- Ms[,i] <- C %*% S_ei %*% fit$y
    if(i %% 1000 ==0) print(range(M))
    # Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fit$y))
    Ss[,((i-1)*ne+1):(i*ne)] <- t(mvtnorm::rmvnorm(n = ne, mean = M, sigma = C, method = "eigen"))
    # Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
    # Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
  }
  
  m <- rowMeans(Ss) #- y_avg
  ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
  lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
  ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
  lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg
  
  if(exists("seed") & exists("ytrue")) write.csv(logs_sample(y = ytrue, dat = Ss), file = paste0("csv/logS/notap/logscore_",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  if(exists("y_hi")) emp_cover <- round(mean(y_hi > lb & y_hi < ub),3)
  if(exists("y_hi")) emp_cover99 <- round(mean(y_hi > lbb & y_hi < ubb),3)
  
  if(exists("ytrue")) emp_cover <- round(mean(ytrue > lb & ytrue < ub),3)
  if(exists("ytrue")) emp_cover99 <- round(mean(ytrue > lbb & ytrue < ubb),3)
  if(exists("seed") & exists("ytrue")) write.csv(mean((m-ytrue)^2), file = paste0("csv/MSE/notap/",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  if(exists("seed") & exists("ytrue")) write.csv(c(emp_cover, emp_cover99), file = paste0("csv/emp_cover_",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  plot(fit$x, fit$y, type="n",
       ylim = range(c(m, lb, ub, lbb, ubb, Y)), 
       main = paste0("est both,", if(exists("emp_cover")){emp_cover}, " ", if(exists("emp_cover99")){emp_cover99}))
  
  for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  if(exists("y_hi")) lines(fit$x, y_hi, lwd=1.5, col="red")
  if(exists("ytrue")) lines(fit$x, ytrue, lwd=1.5, col="red")
  lines(fit$x, fit$y, lwd=1.5, lty=2)
  lines(fit$x, m , col="blue")
  lines(fit$x, lb , col="blue", lty=2)
  lines(fit$x, ub , col="blue", lty=2)
  lines(fit$x, lbb , col="darkblue", lty=2)
  lines(fit$x, ubb , col="darkblue", lty=2)
  legend("bottomright", legend = c("true", "sample avg", "UQ"), 
         col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
  # Using law of total expectation and variance instead
  mu_bar <- rowMeans(Ms)
  C_bar <- matrix(rowMeans(Cs),  length(fit$x),  length(fit$x))
  cov_mut <- cov(t(Ms))
  
  plot(x, mu_bar, type="n", main = "law of total E & V; MPI version", col="blue")
  # for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  lines(x, mu_bar, col="blue")
  if(exists("y_hi")) lines(x, y_hi, col="red")
  if(exists("ytrue")) lines(x, ytrue, col="red")
  lines(x, fit$y, lty=2)
  lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  legend("bottomleft", legend = c("true", "sample avg", "UQ"), 
         col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
}

# Plot simulations of truth given the data, differing tau
plot.true.tau <- function(fit, S_e = fit$Sigma_hat, tolpower = -10, 
                          Y = parent.frame()$Y, unif_tau = 0.1, nrun = nrow(Y)){

  Cs <- matrix(NA, length(fit$x)^2, fit$nmcmc)
  Ss <- St <- Sw <- Sx <- Ms <- matrix(NA, length(fit$x), fit$nmcmc)
  
  # move the next four lines into the for loop if you want to do draws for tau
  print("this fn currently assumes a constant tau")
  
  S_e <- S_e * bohman(plgp:::distance(parent.frame()$x), unif_tau)
  # S_e <- S_e * bohman(plgp:::distance(parent.frame()$x), ifelse(unif_tau, runif(1), exp(rnorm(1))))
  S_ei <- matrix.Moore.Penrose2(S_e, tolp = tolpower)
  S_ei <- (S_ei+t(S_ei))/2
  for (i in 1:fit$nmcmc){
      
     if( i %% 2500 == 0) print(i)
    theta <- fit$theta_y[i]
    # tau2 <- fit$tau2[i]
    # t2S_ei <- 1/tau2 * S_ei
    g <- fit$g[i]
    W <- fit$w[[i]]
    dw <- sq_dist(W)
    if(v==999){
      S_si <- solve(Exp2Fun(distmat = dw, covparms = c(1, theta, g))) #1/tau2 * 
    } else {
      S_si <- solve(MaternFun(distmat = dw, covparms = c(1, theta, g, v))) #1/tau2 * 
    }
    S_si <- (S_si+t(S_si))/2
    C <- Cs[,i] <- solve(S_si + S_ei)
    C <- (C+t(C))/2
    M <- Ms[,i] <- C %*% S_ei %*% fit$y
    if(i %% 1000 ==0) print(range(M))
    Ss[,i] <- M + matrix.sqrt(C)%*%rnorm(length(fit$y))
    # St[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "eigen")
    # Sw[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "svd")
    # Sx[,i] <- mvtnorm::rmvnorm(n = 1, mean = M, sigma = C, method = "chol")
  }
  
  if(exists("seed") & exists("ytrue")) write.csv(logs_sample(y = ytrue, dat = Ss), file = paste0("csv/logS/logscore_",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  m <- rowMeans(Ss) #- y_avg
  ub <- apply(Ss, 1, function(x){quantile(x,0.975)}) #- y_avg
  lb <- apply(Ss, 1, function(x){quantile(x,0.025)}) #- y_avg
  ubb <- apply(Ss, 1, function(x){quantile(x,0.995)}) #- y_avg
  lbb <- apply(Ss, 1, function(x){quantile(x,0.005)}) #- y_avg
  
  if(exists("y_hi")) emp_cover <- round(mean(y_hi > lb & y_hi < ub),3)
  if(exists("y_hi")) emp_cover99 <- round(mean(y_hi > lbb & y_hi < ubb),3)
  
  if(exists("ytrue")) emp_cover <- round(mean(ytrue > lb & ytrue < ub),3)
  if(exists("ytrue")) emp_cover99 <- round(mean(ytrue > lbb & ytrue < ubb),3)
  if(exists("seed") & exists("ytrue")) write.csv(mean((m-ytrue)^2), file = paste0("csv/MSE/tap/",one_layer,"_",seed,"_",nmcmc,".csv"))
  
  if(exists("seed") & exists("ytrue")) write.csv(c(emp_cover, emp_cover99), file = paste0("csv/tap/emp_cover_",one_layer,"_",seed,"_",nmcmc,".csv"))

  plot(fit$x, fit$y, type="n",
       ylim = range(c(m, lb, ub, lbb, ubb, Y)), 
       main = paste0("est both,",paste0("taper \n",unif_tau,"=tau\n"), 
                     if(exists("emp_cover")){emp_cover}, " ", if(exists("emp_cover99")){emp_cover99}))
  
  for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  if(exists("y_hi")) lines(fit$x, y_hi, lwd=1.5, col="red")
  if(exists("ytrue")) lines(fit$x, ytrue, lwd=1.5, col="red")
  lines(fit$x, fit$y, lwd=1.5, lty=2)
  lines(fit$x, m , col="blue")
  lines(fit$x, lb , col="blue", lty=2)
  lines(fit$x, ub , col="blue", lty=2)
  lines(fit$x, lbb , col="darkblue", lty=2)
  lines(fit$x, ubb , col="darkblue", lty=2)
  if(exists("ytrue"))  legend("bottomright", legend = c("true", "sample avg", "UQ"), 
                              col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  if(exists("y_hi"))  legend("bottomright", legend = c("hi_res", "sample avg", "UQ"), 
                             col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
  # Using law of total expectation and variance instead
  mu_bar <- rowMeans(Ms)
  C_bar <- matrix(rowMeans(Cs),  length(fit$x),  length(fit$x))
  cov_mut <- cov(t(Ms))
  
  plot(x, mu_bar, type="n", main = "law of total E & V; MPI version", col="blue",
       ylim=range(Y))
  for (i in 1:nrun) lines(fit$x, Y[i,], col="gray")
  lines(x, mu_bar, col="blue")
  if(exists("y_hi")) lines(x, y_hi, col="red")
  if(exists("ytrue")) lines(x, ytrue, col="red")
  lines(x, fit$y, lty=2)
  lines(x, mu_bar - 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  lines(x, mu_bar + 2*sqrt(diag(C_bar+cov_mut)), col="blue")
  if(exists("ytrue"))  legend("bottomright", legend = c("true", "sample avg", "UQ"), 
                              col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  if(exists("y_hi"))  legend("bottomright", legend = c("hi_res", "sample avg", "UQ"), 
                              col=c("red","black", "blue"), lty=c(1,2,1), lwd=c(2,2,1))
  
}
