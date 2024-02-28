#' Generate values from the posterior predictive distribution
#' This function generate a new dataset with the same dimension of the analyzed one from the posterior predictive distribution
#' @param Z
#' @param Mu
#' @param Sigma
#'
#' @return
#' @export
#'
#' @examples

posterior_predictive_sampling <- function(Z, Mu, Sigma){
  Y.pred <- array(0, dim = c(dim(Z)[1], dim(Mu)[2], dim(Z)[3]))
  for(r in 1:dim(Z)[3]){
    mu.star <- Z[,,r] %*% Mu[,,r]
    if(is.matrix(Sigma)) L.star <- chol(Sigma) else L.star <- chol(Sigma[,,r])
    Y.pred[,,r] <- matrix(rnorm(dim(Z)[1]*dim(Mu)[2]), dim(Z)[1], dim(Mu)[2]) %*% L.star + mu.star
  }
  Y.pred
}
