#' A posteriori log-likelihood computation
#'
#' This function returns the log-likelihood value at each iteration of the MCMC.
#'
#' @param values the output of `perla` function
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#'
#' @return
#' @export
#'
#' @examples
#'

recover.loglikelihood <- function(values, burnin = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  y <- values$y
  Z <- values$Z[,,to.keep]  # n x K x R
  Mu <- values$Mu[,,to.keep]    # K x d x R
  if(is.matrix(values$Sigma)) Sigma <- values$Sigma else Sigma <- Sigma[,,to.keep]
  n <- nrow(y)
  d <- ncol(y)
  val <- numeric(dim(Z)[3])
  pb <- progress_bar$new(total = dim(Z)[3])
  for(r in 1:dim(Z)[3]){
    pb$tick()
    if(is.matrix(Sigma)) S <- Sigma else S <- Sigma[,,r]
    eta <- Z[,,r] %*% Mu[,,r]
    val[r] <- -.5*sum(diag((y - eta) %*% solve(S) %*% t(y - eta)))-
      n*d/2*log(2*pi)-n/2*determinant(S, logarithm = T)$mod
  }
  val
}



model.loglikelihood <- function(y, Z,  # n x K
                                Mu,    # K x d
                                Sigma){
  values <- 0
  n <- nrow(y)
  d <- ncol(y)
  eta <- Z %*% Mu
  -.5*sum(diag((y - eta) %*% solve(Sigma) %*% t(y - eta)))-
    n*d/2*log(2*pi)-n/2*determinant(Sigma, logarithm = T)$mod
}
