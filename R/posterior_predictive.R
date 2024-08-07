#' Generate values from the posterior predictive distribution
#'
#' @description
#' This function generate a new dataset (with the same dimension of the input
#' data) sampling from the posterior predictive distribution.
#'
#' @param values An output of the `perla` function.
#' @param burnin A vector of indexes denoting the MCMC draws to be discarded
#' (default `NULL` means no draws are discarded).
#'
#' @return
#'
#' @export
#'
#' @examples
#'

posterior_predictive_sampling <- function(values, burnin = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  Z <- values$Z[,,to.keep]
  Mu <- values$Mu[,,to.keep]
  Sigma <- values$Sigma[,,to.keep]
  Y.pred <- array(0, dim = c(dim(Z)[1], dim(Mu)[2], dim(Z)[3]))
  for(r in 1:dim(Z)[3]){
    mu.star <- Z[,,r] %*% Mu[,,r]
    L.star <- chol(Sigma[,,r])
    Y.pred[,,r] <- matrix(rnorm(dim(Z)[1]*dim(Mu)[2]), dim(Z)[1], dim(Mu)[2]) %*% L.star + mu.star
  }
  Y.pred
}
