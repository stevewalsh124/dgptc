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
