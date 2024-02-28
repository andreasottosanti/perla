#' Generate values from the posterior predictive distribution
#' This function generate a new dataset with the same dimension of the analyzed one from the posterior predictive distribution
#'
#' @param values the output of `perla` function
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#'
#' @return
#' @export
#'
#' @examples

posterior_predictive_sampling <- function(values, burnin = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  Z <- values$Z[,,to.keep]
  Mu <- values$Mu[,,to.keep]
  if(is.matrix(values$Sigma)) Sigma <- values$Sigma else Sigma <- Sigma[,,to.keep]
  Y.pred <- array(0, dim = c(dim(Z)[1], dim(Mu)[2], dim(Z)[3]))
  for(r in 1:dim(Z)[3]){
    mu.star <- Z[,,r] %*% Mu[,,r]
    if(is.matrix(Sigma)) L.star <- chol(Sigma) else L.star <- chol(Sigma[,,r])
    Y.pred[,,r] <- matrix(rnorm(dim(Z)[1]*dim(Mu)[2]), dim(Z)[1], dim(Mu)[2]) %*% L.star + mu.star
  }
  Y.pred
}
