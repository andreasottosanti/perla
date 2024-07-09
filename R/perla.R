#' PEnalized Regression with Localities Aggregation (PERLA)
#'
#' @description
#' This function estimates a PERLA model using the Markov Chain Monte Carlo
#' algorithm. PERLA is a multivariate Bayesian model that clusters the areas in
#' a territory according to the observed mortality rates of multiple causes of
#' deaths, exploiting also the information of external covariates.
#'
#' @details
#' PERLA incorporates the spatial structure of the data directly into the
#' clustering probabilities by leveraging the stick-breaking formulation of the
#' multinomial distribution. Additionally, it exploits suitable global-local
#' shrinkage priors to ensure that the detection of clusters is driven by
#' concrete differences across mortality levels, while excluding spurious
#' differences. PERLA is estimated with an MCMC algorithm for posterior
#' inference that consists of closed-form Gibbs sampling moves for nearly every
#' model parameter, without requiring complex tuning operations.
#'
#' @import pgdraw
#' @import LaplacesDemon
#' @import sp
#' @import surveillance
#' @import progress
#' @import label.switching
#'
#' @param y Input dataset. It can be either of class `matrix` or
#' `SpatialPolygonsDataFrame`.
#' @param W Adjacency matrix. If not specified, is set to the adjacency matrix
#' extracted from the `SpatialPolygonsDataFrame` data object.
#' @param K Number of clusters.
#' @param R Number of MCMC iterations (default `10^4`).
#' @param prior.rho Type of prior on the rho parameter of the CAR. If `const`,
#' it assumes a fixed value, that is `rho.value`. If `disc`, it assumes that rho
#' can be equal to `0` (with probability `p.spike`) or `rho.value`. If `cont`,
#' it assumes that rho is generated from a mixture of a `dbeta(x,2,18)`
#' (with probability `p.spike`) or a `dbeta(x,18,2)`.
#' @param rho.value Fixed value for the `rho` prior (default `0.99`).
#' @param p.spike
#' @param mean.penalty If `c()`, only the global shrinkage is considered. If
#' `mean.penalty` contains `"c"`, then cluster penalties are added. If
#' `mean.penalty` contains `"d"`, then disease penalties are added. If
#' `mean.penalty` contains `"cd"`, then both cluster and disease penalties are
#' added. If `NULL`, no penalization is considered and the prior of cluster mean
#' vectors are d-variate Gaussian distributions, zero-centered and with
#' covariance matrix `Sigma0`.
#' @param mu0 Vector with values for `mu0` hyperparameter. If NULL (default) the
#' vector components are set equal to 0.
#' @param Sigma0 Matrix with values for `Sigma0` hyperparameter. If NULL
#' (default) the matrix components are set equal to 10.
#' @param burnin A vector of indexes denoting the MCMC draws to be discarded. If
#' `NULL`, then only the starting values are discarded and the algorithm will
#' perform `R+1` iterations.
#' @param initialization It can be either an object of class `perla` or a list
#' of starting points. In the last case, each list element must be named after
#' a parameter to be initialized. Names include `Mu`, `Prob`, `Z`, `Sigma`,
#' `Rho`.
#'
#'
#' @return
#' @export
#'
#' @examples

perla <- function(y, W = NULL, K, R = 10^4,
                  prior.rho = "const",
                  p.spike = .5,
                  rho.value = 0.99,
                  mean.penalty = c(),
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


# Preparation of the objects used to stock the posterior values ----------------

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
  Phi <- Zeta.d <- Zeta.c <- Zeta.cd <- NULL
  if(!is.null(mean.penalty)){
    cat("Applying global penalization;\n")
    Phi <- numeric(R.kept)}
  if("d" %in% mean.penalty){
    cat("Applying feature penalization;\n")
    Zeta.d <- matrix(0, R.kept, d)}
  if("c" %in% mean.penalty){
    cat("Applying cluster penalization;\n")
    Zeta.c <- matrix(0, R.kept, K)}
  if("cd" %in% mean.penalty){
    cat("Applying feature-cluster penalization;\n")
    Zeta.cd <- matrix(0, R.kept, K*d)}
  #cat("              Penalizations\n")
  #cat("----------------------------------------------------")
  #cat("Global   |   Features   |   Clusters   | Features-clusters")


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
  phi.current <- abs(rcauchy(1))
  zeta.d.current <- zeta.c.current <- zeta.cd.current <- NULL
  if("d" %in% mean.penalty) zeta.d.current <- abs(rcauchy(d))
  if("c" %in% mean.penalty) zeta.c.current <- abs(rcauchy(K))
  if("cd" %in% mean.penalty) zeta.cd.current <- abs(rcauchy(K*d))


  detOmega <- determinant((diag(rowSums(W)) - rho.value * W)/tau, logarithm = T)$mod
  acceptance.rho <- numeric(K-1)


  pb <- progress_bar$new(total = R) # progress bar


# Cycle -------------------------------------------------------------------

  for(r in 2:R){
    pb$tick()

    # --update mu
    if(is.null(mean.penalty))
      Mu.current <- update_means(y = y, Z = Z.current, Sigma = Sigma.current, mu0 = mu0, Sigma0 = Sigma0) else {
        results.mu <- update_means_shrink(y = y,
                                          mean.penalty = mean.penalty,
                                          Mu = Mu.current,
                                          Z = Z.current,
                                          Sigma = Sigma.current,
                                          zeta.d = zeta.d.current,
                                          zeta.c = zeta.c.current,
                                          zeta.cd = zeta.cd.current,
                                          phi = phi.current)
        Mu.current <- results.mu$mu
        phi.current <- results.mu$phi
        zeta.d.current <- results.mu$zeta.d
        zeta.c.current <- results.mu$zeta.c
        zeta.cd.current <- results.mu$zeta.cd
      }


    # --update Z
    Z.current <- update_Z(y = y, mu = Mu.current, Sigma = Sigma.current, Psi = Psi.current)

    # --update omega and psi
    updated.psi.omega <- update.psi.omega(psi = Psi.current, omega = Omega.current, Z = Z.current, D = diag(rowSums(W)), W = W, tau = tau, rho = Rho.current, apply.mean.correction = F)
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

    # --update Sigma
    resid <- y - Z.current %*% Mu.current
    Sigma.current <- rinvwishart(nu = d + n, S = diag(d) + t(resid) %*% resid)

    # ---store the values
    if(r %in% tokeep){
      Mu[,,r-max(burnin)] <- Mu.current
      Z[,,r-max(burnin)] <- Z.current
      Prob[,,r-max(burnin)] <- Prob.current
      Rho[r-max(burnin),] <- Rho.current
      Sigma[,,r-max(burnin)] <- Sigma.current
      if(!is.null(mean.penalty)) Phi[r-max(burnin)] <- phi.current
      if("d" %in% mean.penalty) Zeta.d[r-max(burnin),] <- zeta.d.current
      if("c" %in% mean.penalty) Zeta.c[r-max(burnin),] <- zeta.c.current
      if("cd" %in% mean.penalty) Zeta.cd[r-max(burnin),] <- zeta.cd.current
    }
  }

  cat("\n")
  shrinkage.parameters <- list()
  if(!is.null(mean.penalty)) shrinkage.parameters$Phi <- Phi
  if("d" %in% mean.penalty) shrinkage.parameters$Zeta.d <- Zeta.d
  if("c" %in% mean.penalty) shrinkage.parameters$Zeta.c <- Zeta.c
  if("cd" %in% mean.penalty) shrinkage.parameters$Zeta.cd <- Zeta.cd
  results <- list(Mu = Mu, Z = Z, Sigma = Sigma, Prob = Prob, Rho = Rho,
                  shrinkage.parameters = shrinkage.parameters,
                  acceptance.rho = acceptance.rho, y = y)
  class(results) <- "perla"
  return(results)

}
