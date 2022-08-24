matrix.Moore.Penrose <- function(H)
{
  # Computes Moore-Penrose generalized inverse of symmetric matrix using spectral decomposition
  # By M.A.R. Ferreira
  H.eigen = eigen(H)
  inverse.values = rep(0,nrow(H))
  inverse.values[abs(H.eigen$values) > 10^(-10)] = 1/H.eigen$values[abs(H.eigen$values) > 10^(-10)]
  H.MP = H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  H.MP
}

matrix.Moore.Penrose2 <- function(H, tolp = -10)
{
  # Computes Moore-Penrose generalized inverse of symmetric matrix using spectral decomposition
  # By M.A.R. Ferreira
  H.eigen = eigen(H)
  inverse.values = rep(0,nrow(H))
  inverse.values[(H.eigen$values) > 10^(tolp)] = 1/H.eigen$values[(H.eigen$values) > 10^(tolp)]
  H.MP = H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors)
  H.MP
}

matrix.sqrt <- function(H)
{
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-20)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  
  
  H.sqrt
}
