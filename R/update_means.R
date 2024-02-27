update_centroids <- function(y, Z, Sigma, mu0, Sigma0){
  if(is.vector(Z)) Z <- matrix(Z, ncol = 1)
  K <- ncol(Z)
  invSigma <- solve(Sigma)
  invSigma0 <- solve(Sigma0)
  mu <- matrix(0, K, ncol(y))
  for(k in 1:K){
    n_k <- sum(Z[,k] == 1)
    if(n_k > 0){
      if(n_k > 1){
        Var <- solve(invSigma0 + n_k * invSigma)
        Mean <- Var %*% (invSigma0 %*% mu0 +invSigma %*% colSums(y[Z[,k] == 1,]))
        mu[k,] <- t(chol(Var)) %*% rnorm(ncol(y)) + Mean
      } else {
        Var <- solve(invSigma0 + n_k * invSigma)
        Mean <- Var %*% (invSigma0 %*% mu0 +invSigma %*% y[Z[,k] == 1,])
        mu[k,] <- t(chol(Var)) %*% rnorm(ncol(y)) + Mean
      }
    } else {
      mu[k,] <- t(chol(Sigma0)) %*% rnorm(ncol(y)) + mu0
    }
  }
  mu
}
