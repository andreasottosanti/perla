#' Generate samples from the posterior predictive distribution
#'
#' @description
#' This function generates a collection of datasets (each of the same size of the input
#' data) sampling from the posterior predictive distribution.
#'
#' @param x An object of class `perla`.
#' @param burnin A vector of indexes denoting the MCMC draws to be discarded
#' (default `NULL` means no draw to be discarded).
#'
#' @return An array of dimension `n x d x R`, where `n` is the number of areas in the map, `d` is the number of variables, and `R` is the number of MCMC iterations (minus the burnin, is given).
#'
#' @export
#'
#' @examples
#'

posterior_predictive_sampling <- function(x, burnin = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(x$Z)[3] else
    to.keep <- setdiff(1:dim(x$Z)[3], burnin)
  Z <- x$Z[,,to.keep]
  Mu <- x$Mu[,,to.keep]
  Sigma <- x$Sigma[,,to.keep]
  Y.pred <- array(0, dim = c(dim(Z)[1], dim(Mu)[2], dim(Z)[3]))
  for(r in 1:dim(Z)[3]){
    mu.star <- Z[,,r] %*% Mu[,,r]
    L.star <- chol(Sigma[,,r])
    Y.pred[,,r] <- matrix(rnorm(dim(Z)[1]*dim(Mu)[2]), dim(Z)[1], dim(Mu)[2]) %*% L.star + mu.star
  }
  Y.pred
}
