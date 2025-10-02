update_responsabilities <- function(y, mu, Sigma, Psi){
  if(is.vector(Psi)) Psi <- matrix(Psi, ncol = 1)
  n <- nrow(y)
  K <- ncol(Psi)+1
  Prob <- Z <- matrix(0, n, K)
  pp <- convert.to.probabilities(Psi)
  for(k in 1:K){
    q1 <- y %*% solve(Sigma) %*% mu[k,] - as.vector(.5 * t(mu[k,]) %*% solve(Sigma) %*% mu[k,])
    Prob[,k] <- exp(q1)*pp[,k]
  }
  return(Prob)
}

update_Z <- function(responsabilities){
  n <- nrow(responsabilities)
  K <- ncol(responsabilities)
  Z <- matrix(0, n, K)
  for(i in 1:n){
    z <- sample(x = 1:K, size = 1, prob = responsabilities[i,])
    Z[i,z] <- 1
  }
  Z
}


convert.to.probabilities <- function(psi){
  if(is.vector(psi)) psi <- matrix(psi, nrow = 1)
  prob <- matrix(0, nrow(psi), ncol(psi)+1)
  prob[,1] <- exp(psi[,1])/(1+exp(psi[,1]))
  if(ncol(psi) > 1){
    for(k in 2:ncol(psi)){
      prob.tilde <- exp(psi[,k])/(1+exp(psi[,k]))
      if(k == 2)
        prob[,2] <- (1-prob[,1])*prob.tilde else
          prob[,k] <- (1-rowSums(prob[,1:(k-1)]))*prob.tilde
    }
    prob[,ncol(psi)+1] <- (1-rowSums(prob[,1:ncol(psi)]))
  } else {
    prob[,ncol(psi)+1] <- 1-prob[,1]
  }
  prob
}

convert.to.stickbreaking <- function(p){
  K <- length(p)
  prob <- numeric(K)
  for(k in 1:K)
    if(k == 1) prob[k] <- p[k] else
      prob[k] <- p[k]/(1-sum(p[1:(k-1)]))
  prob
}



