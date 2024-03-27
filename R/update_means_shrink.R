update_means_shrink <- function(y, mean.penalty, Z, Mu, Sigma, zeta.d, zeta.c, zeta.cd, phi){
  if(is.vector(Z)) Z <- matrix(Z, ncol = 1)
  K <- ncol(Z)
  d <- ncol(y)
  invSigma <- solve(Sigma)
  mu <- matrix(0, K, ncol(y))
  if(is.null(zeta.d)) zeta.d.current <- rep(1, d) else zeta.d.current <- zeta.d
  if(is.null(zeta.c)) zeta.c.current <- rep(1, K) else zeta.c.current <- zeta.c
  if(is.null(zeta.cd)) zeta.cd.current <- rep(1, K*d) else zeta.cd.current <- zeta.cd
  zeta.cd.current <- matrix(zeta.cd.current, K, d)
  phi.current <- phi
  #if(include.Sigma){
  #  L <- chol(Sigma)
  #  Mu <- Mu %*% solve(L)
  #}

# update zetas ------------------------------------------------------------
  if("d" %in% mean.penalty){
    for(j in 1:d){
      a_zeta_j <- rexp(1, rate = zeta.d.current[j]+1)
      zeta.d.current[j] <- rgamma(1, (K+1)/2, phi.current*sum(zeta.cd.current[,j] * zeta.c.current * Mu[,j]^2)/2 + a_zeta_j)
    }
  }
  if("c" %in% mean.penalty | "cd" %in% mean.penalty){
    for(k in 1:K){
      if("c" %in% mean.penalty){
        a_zeta_k <- rexp(1, rate = zeta.c.current[k]+1)
        zeta.c.current[k] <- rgamma(1, (d+1)/2, phi.current*sum(zeta.cd.current[k,] * zeta.d.current * Mu[k,]^2)/2 + a_zeta_k)
      }
      if("cd" %in% mean.penalty){
        for(j in 1:d){
          a_zeta_jk <- rexp(1, rate = zeta.cd.current[k,j]+1)
          zeta.cd.current[k,j] <- rgamma(1, 1, phi.current * zeta.d.current[j] * zeta.c.current[k] * Mu[k,j]^2/2 + a_zeta_jk)
        }
      }
    }
  }

# update phi --------------------------------------------------------------
  a_phi <- rexp(1, rate = phi.current+1)
  phi.current <- rgamma(1, (K*d+1)/2, (t(zeta.c.current) %*% (Mu^2 * zeta.cd.current) %*% zeta.d.current)/2 + a_phi)
  #phi.current <- rgamma(1, (K*d+1)/2, sum(zeta.current * (t(Mu^2) %*% rep(1,K)))/2 + a_phi)

# update mu ---------------------------------------------------------------
  for(k in 1:K){
    n_k <- sum(Z[,k] == 1)
    if(n_k > 0){
      if(n_k > 1){
        Var <- solve(phi.current * zeta.c.current[k] * diag(zeta.d.current * zeta.cd.current[k,]) + n_k * invSigma)
        Mean <- Var %*% (invSigma %*% colSums(y[Z[,k] == 1,]))
        mu[k,] <- t(chol(Var)) %*% rnorm(ncol(y)) + Mean
      } else {
        Var <- solve(phi.current * zeta.c.current[k] * diag(zeta.d.current * zeta.cd.current[k,]) + n_k * invSigma)
        Mean <- Var %*% (invSigma %*% y[Z[,k] == 1,])
        mu[k,] <- t(chol(Var)) %*% rnorm(ncol(y)) + Mean
      }
    } else {
      mu[k,] <- rnorm(d, sd = sqrt(1/(phi.current * zeta.c.current[k] * zeta.d.current * zeta.cd.current[k,])))
    }
  }
  return(list(mu = mu, zeta.c = zeta.c.current, zeta.d = zeta.d.current, zeta.cd = as.vector(zeta.cd.current), phi = phi.current))
}


