#' Posterior update of the CAR variance parameters via Metropolis-Hastings
#'
#' This function update the CAR variance parameter drawing proposing values from a log-Normal distribution
#'
#' @importFrom emulator quad.diag
#'
#' @param tau_current the current parameter value.
#' @param Psi a matrix of dimension `n x (K-1)`, where `n` is the number of spatial areas in the analysed map, and `K` is the number of clusters assumed.
#' @param D the diagonal matrix containing the number of neighbours.
#' @param W the proximity matrix.
#' @param Rho a vector of length `K-1` containing the CAR \eqn{(\rho_k)_k} values.
#' @param shape shape value of the inverse gamma assumed a priori.
#' @param scale scale value of the inverse gamma assumed a priori.
#' @param eps a tuning parameter (default `1`). Proposed values are generated from `rlnorm` with `meanlog = log(tau_current)` and `sdlog = eps`.
#'
#' @return A list of two elements:
#'
#' - `value` contains the value simulated by the algorithm;
#' - `accepted` is a value that equals 1 if the returned value is different from `tau_current`, and 0 otherwise.
#'
update_tau_metropolis <- function(tau_current, Psi, D, W, Rho, shape, scale, eps = 1){
  shape_post <- ncol(Psi)*nrow(Psi)/2 + shape
  output_values <- numeric(ncol(Psi))
  quadr.form <- quad.diag(M = D, x = Psi) - Rho * quad.diag(M = W, x = Psi)
  scale_post <- sum(quadr.form)/2 + scale
  log_tau_star <- rnorm(1, log(tau_current), eps)
  accepted <- 0
  A <- exp(dinvgamma(exp(log_tau_star), shape = shape_post, scale = scale_post, log = T) -
    dinvgamma(tau_current, shape = shape_post, scale = scale_post, log = T) +
    dlnorm(tau_current, log_tau_star, eps, log = T) -
    dlnorm(exp(log_tau_star), log(tau_current), eps, log = T)
  )
  if(runif(1) <= A){
    tau_current <- exp(log_tau_star)
    accepted <- 1
  }
  return(list(value = tau_current, accepted = accepted))
}
