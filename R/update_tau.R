#' Posterior update of the CAR variance parameters via conjugate priors
#'
#' This function update the CAR variance parameter drawing from its conditional posterior distributions, which is Inverge Gamma.
#'
#' @importFrom emulator quad.diag
#'
#' @param Psi a matrix of dimension `n x (K-1)`, where `n` is the number of spatial areas in the analysed map, and `K` is the number of clusters assumed.
#' @param D the diagonal matrix containing the number of neighbours.
#' @param W the proximity matrix.
#' @param Rho a vector of length `K-1` containing the CAR \eqn{(\rho_k)_k} values.
#' @param shape shape value of the inverse gamma assumed a priori.
#' @param scale scale value of the inverse gamma assumed a priori.
#'
#' @return a value generated from the conditional posterior distribution of \eqn{\tau}.
#'

update_tau <- function(Psi, D, W, Rho, shape, scale){
  shape_post <- ncol(Psi)*nrow(Psi)/2 + shape
  output_values <- numeric(ncol(Psi))
  quadr.form <- quad.diag(M = D, x = Psi) - Rho * quad.diag(M = W, x = Psi)
  scale_post <- sum(quadr.form)/2 + scale
  return(rinvgamma(n = 1, shape = shape_post, scale = scale_post))
}
