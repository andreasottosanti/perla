library(LaplacesDemon)

logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))

proposal.rho <- function(rho, mean = 0, sd = 1){
  1/(rho*(1-rho))*dnorm(logit(rho), mean, sd)*I(rho > 0 & rho < 1)
}

prior.rho <- function(rho,
                      lambda_spike = 3,
                      lambda_slab = 7,
                      a_spike = .5,
                      b_spike = 18,
                      a_slab = 18,
                      b_slab = 2){
  lambda_spike/(lambda_spike+lambda_slab)*dbeta(rho, a_spike, b_spike) + lambda_slab/(lambda_spike+lambda_slab)*dbeta(rho, a_slab, b_slab)
}

log.posterior.rho <- function(rho, psi, tau = 1, W, lambda_spike = 3,
                              lambda_slab = 7,
                              a_spike = .5,
                              b_spike = 18,
                              a_slab = 18,
                              b_slab = 2){
  D <- rowSums(W)
  mean.vector <- rho * (W %*% psi) / D
  sd.vector <- sqrt(tau / D)
  sum(dnorm(psi, mean = mean.vector, sd = sd.vector, log = T)) +
    log(prior.rho(rho, lambda_spike, lambda_slab, a_spike, b_spike, a_slab, b_slab))
}


metropolis.rho <- function(rho.old, psi, tau = 1, W,
                           lambda_spike = 3,
                           lambda_slab = 7,
                           a_spike = .5,
                           b_spike = 18,
                           a_slab = 18,
                           b_slab = 2,
                           sigma.rho = 1){
  r.star <- rnorm(1, logit(rho.old), sigma.rho)
  rho.star <- expit(r.star)
  proposal.num <- proposal.rho(rho = rho.old, mean = logit(rho.star), sd = sigma.rho)
  proposal.den <- proposal.rho(rho = rho.star, mean = logit(rho.old), sd = sigma.rho)
  accepted <- 0
  A <- exp(log.posterior.rho(rho = rho.star, psi = psi, tau = tau, W = W, lambda_spike, lambda_slab, a_spike, b_spike, a_slab, b_slab)-
             log.posterior.rho(rho = rho.old, psi = psi, tau = tau, W = W, lambda_spike, lambda_slab, a_spike, b_spike, a_slab, b_slab))*
    proposal.num/proposal.den
  if(runif(1) <= A){
    accepted <- 1
    rho.old <- rho.star
  }
  return(list(rho = rho.old, accepted = accepted))
}

#' Metropolis-Hastrings update of the CAR parameter
#'
#' This function receives a value of rho and generates a new value from its posterior distribution given by the spike & slab prior.
#' It uses a random walk Metropolis-Hastrings move.
#'
#' @import LaplacesDemon
#' @export
#'
#' @param rho.old
#' @param psi
#' @param tau
#' @param W
#' @param lambda_spike
#' @param lambda_slab
#' @param a_spike
#' @param b_spike
#' @param a_slab
#' @param b_slab
#' @param sigma.rho
#'
#' @return
#'
#' @examples
#'


update_rho_continuous <- function(rho.old, # a vector of length K-1
                       psi,  # the n x (K-1) matrix
                       tau = 1, W,
                       lambda_spike = 3,
                       lambda_slab = 7,
                       a_spike = .5,
                       b_spike = 18,
                       a_slab = 18,
                       b_slab = 2,
                       sigma.rho = 1){
  if(length(sigma.rho) == 1 & !is.vector(psi)) sigma.rho <- rep(sigma.rho, ncol(psi))
  if(is.vector(psi)) psi <- matrix(psi, ncol = 1)
  acceptances <- rhos <- numeric(ncol(psi))
  for(k in 1:ncol(psi)){
    rho.metropolis.results <- metropolis.rho(rho.old = rho.old[k], psi = psi[,k], tau = tau, W = W, lambda_spike, lambda_slab, a_spike, b_spike, a_slab, b_slab, sigma.rho = sigma.rho[k])
    acceptances[k] <- acceptances[k]+rho.metropolis.results$accepted
    rhos[k] <- rho.metropolis.results$rho
  }
  return(list(rho = rhos, accepted = acceptances))
}


