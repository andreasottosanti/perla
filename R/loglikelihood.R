library(LaplacesDemon)

model.loglikelihood <- function(y, Z,  # n x K
                                Mu,    # K x d
                                Sigma){
  value <- 0
  n <- nrow(y)
  d <- ncol(y)
  eta <- Z %*% Mu
  -.5*sum(diag((y - eta) %*% solve(Sigma) %*% t(y - eta)))-
    n*d/2*log(2*pi)-n/2*determinant(Sigma, logarithm = T)$mod
}

recover.loglikelihood <- function(y,
                                  Z,  # n x K x R
                                  Mu,    # K x d x R
                                  Sigma){
  value <- 0
  n <- nrow(y)
  d <- ncol(y)
  values <- numeric(dim(Z)[3])
  cat("|")
  for(r in 1:dim(Z)[3]){
    if(r %% 1000 == 0) cat("=")
    if(is.matrix(Sigma)) S <- Sigma else S <- Sigma[,,r]
    eta <- Z[,,r] %*% Mu[,,r]
    values[r] <- -.5*sum(diag((y - eta) %*% solve(S) %*% t(y - eta)))-
      n*d/2*log(2*pi)-n/2*determinant(S, logarithm = T)$mod
  }
  cat("|\n")
  values
}
