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


