#' PEnalized Regression with Local Aggregation (PERLA)
#'
#' This function estimates the PERLA model using the Markov Chain Monte Carlo algorithm.
#'
#' @import pgdraw
#' @import LaplacesDemon
#' @import sp
#' @import surveillance
#' @import progress
#' @import label.switching
#'
#' @param y the input dataset. It can be either of class `matrix` or `SpatialPolygonsDataFrame`
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
#' @param initialization it can be either an object of class `perla` or a list with starting point. In the last case, each list element must be named after the parameter to be initialized. Names include `Mu`, `Prob`, `Z`, `Sigma`, `Rho`)
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded. If `NULL`, then only the starting values are discarded and the algorithm will perform `R+1` iterations.
#'
#' @return
#' @export
#'
#' @examples

perla <- function(y, W = NULL, K, R = 10^4,
                  prior.rho = "const",
                  p.spike = .5,
                  rho.value = 0.99,
                  mu0 = NULL,
                  Sigma0 = NULL,
                  burnin = NULL,
                  initialization = NULL){
  if("SpatialPolygonsDataFrame" %in% class(y)){
    if(is.null(W)) W <- poly2adjmat(y)
    y <- as.matrix(y@data)
  }
  n <- nrow(y)
  d <- ncol(y)
  tau <- 1
  if(p.spike < 0 | p.spike > 1) stop("p.spike is not a probability")

# Preparation of the objects used to stock the draws ----------------------

  if(is.null(burnin)){
    burnin <- 1
    R <- R+1
  }
  tokeep <- setdiff(1:R, burnin)
  R.kept <- length(tokeep)
  Rho <- matrix(0, R.kept, K-1)
  Mu <- array(0, dim = c(K, d, R.kept))
  Z <- Prob <- array(0, dim = c(n, K, R.kept))
  Sigma <- array(0, dim = c(d, d, R.kept))

# Default hyperparameters -------------------------------------------------
  if(is.null(mu0)) mu0 <- rep(0, d)
  if(is.null(Sigma0)){
    if(d > 1) Sigma0 <- diag(rep(10, d)) else
      Sigma0 <- 10}

# Inizialization ----------------------------------------------------------
  if(!is.null(initialization) & class(initialization) == "perla"){
    max.iter.inizialization <- dim(initialization$Z)[3]
    Z.current <- initialization$Z[,,max.iter.inizialization]
    Rho.current <- initialization$Rho[max.iter.inizialization,]
    Mu.current <- initialization$Mu[,,max.iter.inizialization]
    Sigma.current <- initialization$Sigma[,,max.iter.inizialization]
  } else {
    if(is.null(initialization$Z)){
      Z.current <- matrix(0, n, K)
      for(i in 1:n) Z.current[i,sample(1:K, 1)] <- 1
    } else Z.current <- initialization$Z
    if(is.null(initialization$Mu))
      Mu.current <- solve(t(Z.current) %*% Z.current) %*% t(Z.current) %*% y else
        Mu.current <- initialization$Mu
    if(is.null(initialization$Rho))
      Rho.current <- rep(rho.value, K-1) else
        Rho.current <- initialization$Rho
    if(is.null(initialization$Sigma))
      Sigma.current <- cov(y) else
        Sigma.current <- initialization$Sigma
  }
  Psi.current <- matrix(rnorm(n*(K-1)), n, K-1)
  Omega.current <- matrix(rgamma(n*(K-1), 1, 1), n, K-1)

  detOmega <- determinant((diag(rowSums(W)) - rho.value * W)/tau, logarithm = T)$mod
  acceptance.rho <- numeric(K-1)


  pb <- progress_bar$new(total = R) # progress bar

  for(r in 2:R){
    pb$tick()

    # --update mu
    Mu.current <- update_means(y = y, Z = Z.current, Sigma = Sigma.current, mu0 = mu0, Sigma0 = Sigma0)

    # --update Z
    Z.current <- update_Z(y = y, mu = Mu.current, Sigma = Sigma.current, Psi = Psi.current)

    # --update omega and psi
    updated.psi.omega <- update.psi.omega(psi = Psi.current, omega = Omega.current, Z = Z.current, D = D, W = W, tau = tau, rho = Rho.current, apply.mean.correction = F)
    Psi.current <- updated.psi.omega$psi
    Omega.current <- updated.psi.omega$omega
    if(r %in% tokeep) Prob.current <- convert.to.probabilities(Psi.current)

    # --update rho

    # discrete S&S ------------------------------------------------------------
    #rho.updated <- update_rho(rho.old = Rho.current, psi = Psi.current, tau = tau, W = W, lambda1 = p.spike, lambda2 = 1-p.spike,
    #                          p.nonzero = p.nonzero, sigma.rho = sigma.rho)
    #acceptance.rho <- acceptance.rho + rho.updated$accepted
    #Rho.current <- rho.updated$rho

    if(prior.rho == "disc")
      Rho.current <- update_rho_discrete(psi = Psi.current,
                                         tau = tau,
                                         W = W,
                                         lambda.spike = p.spike,
                                         lambda.slab = 1-p.spike,
                                         rho.value = rho.value,
                                         detOmega = detOmega)
    if(prior.rho == "cont"){
      rho.updated <- update_rho_continuous(rho.old = Rho.current,
                                           psi = Psi.current,
                                           tau = tau,
                                           W = W,
                                           lambda_spike = p.spike,
                                           lambda_slab = 1-p.spike, sigma.rho = 3)
      acceptance.rho <- acceptance.rho + rho.updated$accepted
      Rho.current <- rho.updated$rho
    }

    # ---store the values
    if(r %in% tokeep){
      Mu[,,r-max(burnin)] <- Mu.current
      Z[,,r-max(burnin)] <- Z.current
      Prob[,,r-max(burnin)] <- Prob.current
      Rho[r-max(burnin),] <- Rho.current
      Sigma[,,r-max(burnin)] <- Sigma.current
    }
  }

  results <- list(Mu = Mu, Z = Z, Sigma = Sigma, Prob = Prob, Rho = Rho, acceptance.rho = acceptance.rho, y = y)
  class(results) <- "perla"
  return(results)

}
