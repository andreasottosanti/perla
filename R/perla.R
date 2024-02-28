#' PEnalized Regression with Local Aggregation (PERLA)
#'
#' This function estimates the PERLA model using the Markov Chain Monte Carlo algorithm.
#'
#' @import pgdraw
#' @import LaplacesDemon
#' @import sp
#' @import surveillance
#' @import progress
#'
#' @param x can be either a `SpatialPolygonsDataFrame` object or an `n` x `d` matrix, where `n` is the number of spatial areas and `d` is the number of diseases.
#' @param K the number of clusters
#' @param R the number of MCMC iterations (default `10^4`)
#' @param prior.rho set the type of prior on the rho parameter of the CAR. If `const`, it assumes a fixed value, that is `rho.value`.
#' If `disc`, it assumes that rho can be equal to `0` (with probability `p.spike`) or `rho.value`.
#' If `cont`, it assumes that rho is generated from a mixture of a `dbeta(x,2,18)` (with probability `p.spike`) or a `dbeta(x,18,2)`.
#' @param rho.value default `0.99`
#' @param p.spike
#' @param mu0
#' @param Sigma0
#' @param W
#' @param initialization the output of a `perla` function that is used as a starting point of the algorithm
#'
#' @return
#' @export
#'
#' @examples

perla <- function(x, W = NULL, K, R = 10^4,
                  prior.rho = "const",
                  p.spike = .5,
                  rho.value = 0.99,
                  mu0 = NULL,
                  Sigma0 = NULL,
                  initialization = NULL){
  if("SpatialPolygonsDataFrame" %in% class(x)){
    if(is.null(W)) W <- poly2adjmat(map)
    x <- x@data
  }
  n <- nrow(y)
  d <- ncol(y)
  tau <- 1
  if(p.spike < 0 | p.spike > 1) stop("p.spike is not a probability")
  lambda.spike <- p.spike
  lambda.slab <- 1-p.spike

  Rho <- matrix(0, R, K-1)
  Psi <- Omega <- array(0, dim = c(n, K-1, R))
  Mu <- array(0, dim = c(K, d, R))
  Z <- Prob <- array(0, dim = c(n, K, R))

# Default hyperparameters -------------------------------------------------
  if(is.null(mu0)) mu0 <- rep(0, d)
  if(is.null(Sigma0)){
    if(d > 1) Sigma0 <- diag(rep(10, d)) else
      Sigma0 <- 10}

# Inizialization ----------------------------------------------------------
  if(is.null(initialization)){
    for(i in 1:n){
      v <- sample(1:K, 1)
      Z[i,v,1] <- 1
    }
    Rho[1,] <- rho.value
  } else {
    max.iter.inizialization <- dim(initialization$Z)[3]
    Z[,,1] <- initialization$Z[,,max.iter.inizialization]
    Rho[1,] <- initialization$Rho[max.iter.inizialization,]
    Mu[,,1] <- initialization$Mu[,,max.iter.inizialization]
  }
  Psi[,,1] <- matrix(rnorm(n*(K-1)), n, K-1)
  Omega[,,1] <- rgamma(n*(K-1), 1, 1)

  detOmega <- determinant((diag(rowSums(W)) - rho.value * W)/tau, logarithm = T)$mod
  acceptance.rho <- numeric(K-1)

  pb <- progress::progress_bar$new(total = R) # progress bar
  for(r in 2:R){
    progress::pb$tick()

    # ---overwrite
    Mu[,,r] <- Mu[,,r-1]
    Z[,,r] <- Z[,,r-1]
    Psi[,,r] <- Psi[,,r-1]
    Omega[,,r] <- Omega[,,r-1]
    Rho[r,] <- Rho[r-1,]

    # --update mu
    Mu[,,r] <- update_means(y = y, Z = Z[,,r], Sigma = Sigma, mu0 = mu0, Sigma0 = Sigma0)

    # --update Z
    Z[,,r] <- update_Z(y = y, mu = Mu[,,r], Sigma = Sigma, Psi = Psi[,,r])

    # --update omega and psi
    updated.psi.omega <- update.psi.omega(psi = Psi[,,r], omega = Omega[,,r], Z = Z[,,r], D = D, W = W, tau = tau, rho = Rho[r,], apply.mean.correction = F)
    Psi[,,r] <- updated.psi.omega$psi
    Omega[,,r] <- updated.psi.omega$omega
    Prob[,,r] <- convert.to.probabilities(Psi[,,r])

    # --update rho

    # discrete S&S ------------------------------------------------------------
    #rho.updated <- update_rho(rho.old = Rho[r,], psi = Psi[,,r], tau = tau, W = W, lambda1 = lambda.spike, lambda2 = lambda.slab,
    #                          p.nonzero = p.nonzero, sigma.rho = sigma.rho)
    #acceptance.rho <- acceptance.rho + rho.updated$accepted
    #Rho[r,] <- rho.updated$rho


    # rho ---------------------------------------------------------
    if(prior.rho == "disc")
      Rho[r,] <- update_rho_discrete(psi = Psi[,,r], tau = tau, W = W, lambda.spike = lambda.spike, lambda.slab = lambda.slab, rho.value = rho.value, detOmega = detOmega)
    if(prior.rho == "cont"){
      rho.updated <- update_rho_continuous(rho.old = Rho[r,], psi = Psi[,,r], tau = tau, W = W,
                                           lambda_spike = lambda.spike,
                                           lambda_slab = lambda.slab, sigma.rho = 3)
      acceptance.rho <- acceptance.rho + rho.updated$accepted
      Rho[r,] <- rho.updated$rho
    }
  }

  return(list(Mu = Mu, Z = Z, Sigma = Sigma, Prob = Prob, Rho = Rho, acceptance.rho = acceptance.rho, y = y))

}
