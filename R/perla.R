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
#' @param W Adjacency matrix. It is required if `y` is of class `matrix`. If `y` is a `SpatialPolygonsDataFrame` object and `W` is not specified, then it is set equal to the adjacency matrix
#' extracted from `y`.
#' @param K Number of clusters.
#' @param R Number of MCMC iterations performed (default `10^4`).
#' @param X a `data.frame` object containing the variables to use as covariates.
#' @param prior.rho Type of prior on the `rho` parameter of the CAR. If the
#' parameter is set to `const`, it is assumed as fixed value (`rho.value`). If
#' the parameter is set to `disc`, `rho` is allowed to take value `0` (with probability `p.spike`) or `rho.value` (with probability `1 - p.spike`).
#' If `cont`, it assumes that rho is generated from a `dbeta(x, 2, 18)` with probability `p.spike` or from a `dbeta(x, 18, 2)` with probability `1 - p.spike`.
#' @param rho.value Fixed value for the `rho` prior (default `0.99`).
#' @param p.spike If the parameter `prior.rho`  is set to `disc`, `p.spike`
#' indicates the probability of `rho` being `0`.
#' @param mean.penalty A vector specifying the shrinkage factors for the elements of `Mu`.
#' It can include the following elements:
#' - `1`: applies global shrinkage;
#' - `"c"`: applies cluster-specific shrinkage;
#' - `"d"`: applies feature-specific shrinkage;
#' - `"cd"`: applies cluster-feature-specific shrinkage.
#'
#' If `"c"`, `"d"`, or `"cd"` are used, it is not necessary to include `1` in `mean.penalty`.
#' Note that using only `"cd"` corresponds to applying the horseshoe prior.
#' If `NULL`, no penalization is applied, and the prior for the cluster mean vectors is a
#' d-variate Gaussian distribution centered at `mu0` with covariance matrix `Sigma0`.
#' @param mu0 Vector with values for `mu0` hyperparameter. If NULL (default) the
#' vector components are set equal to 0.
#' @param Sigma0 Matrix with values for `Sigma0` hyperparameter. If NULL
#' (default) the matrix components are set equal to 10.
#' @param burnin A vector of indices specifying the MCMC draws to discard from the `R` iterations.
#' The resulting chains will have a length of `R - length(burnin)`.
#' If set to `NULL`, only the initial values are discarded, and the algorithm
#' will run for `R + 1` iterations.
#' @param initialization It can be either an object of class `perla` or a list
#' of starting points. In the last case, each list element must be named after
#' the parameters to be initialized. Names include `Mu`, `Prob`, `Z`, `Sigma`, `Rho`.
#' @param prior.tau Type of prior on the `tau` parameter of the CAR. If the
#' parameter is set to `const`, it is assumed as fixed value (`tau.value`). If
#' the parameter is set to `conj` or `metropolis`, then the parameter is updated via MCMC, considering the conjugate prior `dinvgamma(x, tau.hyperpar[1], tau.hyperpar[2])`.
#' In the first case, a Gibbs sampling is used, otherwise the sampling is performed via Metropolis-Hastings algorithm, generating candidate values from a log-Normal distribution.
#' @param tau.value A value for the marginal variance parameter of the CAR model (default is 1).
#' If `prior.tau == 'const'`, then `tau.value` remains fixed, while if `prior.tau == 'conj'`, then `tau.value` is used to initialise the MCMC algorithm.
#' @param tau.hyperpar A vector of length 2 containing the shape and the scale of the inverse gamma distribution used as a prior for tau.
#' If `tau.value == 'const'`, this is ignored.
#' @param tau.tuning A tuning parameter used when `tau.hyperpar == metropolis'` (default `1`). Proposed values are generated from `rlnorm` with `meanlog = log(tau_current)` and `sdlog = tau.tuning`.
#' If `tau.value != 'metropolis'`, this is ignored.
#' @return
#' An object of class `perla`.
#'
#' @details
#' The model parameters and the clustering labels are initialised at random.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # running Perla without shrinkage parameters
#' perla(y = west_states_data, W = west_states_W, K = 3)
#'
#' # running Perla with only the global shrinkage parameter
#' perla(y = west_states_data, W = west_states_W, K = 3, mean.penalty = c(1))
#'
#' # running Perla with the global shrinkage parameter and the cluster-specific shrinkage parameters
#' perla(y = west_states_data, W = west_states_W, K = 3, mean.penalty = c("c"))
#' }
#'

perla <- function(y, W = NULL, K, R = 10^4,
                  X = NULL,
                  prior.rho = "const",
                  rho.value = 0.99,
                  p.spike = .5,
                  prior.tau = "const",
                  tau.value = 1,
                  tau.hyperpar = c(2.1, 3.1),
                  tau.tuning = 1,
                  mean.penalty = c(1),
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


# some checks -------------------------------------------------------------
  if(p.spike < 0 | p.spike > 1) stop("p.spike is not a probability")
  if(!(prior.rho %in% c("const", "disc", "cont"))) stop("prior.rho can be 'const', 'disc', or 'cont'")
  if(!(prior.tau %in% c("const", "conj", "metropolis"))) stop("prior.tau can be 'const', 'conj', or 'metropolis'")
  if(length(tau.hyperpar) != 2) stop("tau.hyperpar must be a vector of length 2")


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
  Tau <- numeric(R.kept)
  Phi <- Zeta.d <- Zeta.c <- Zeta.cd <- NULL


# display the shrinkage factors selected ----------------------------------

  unicode_global <- unicode_features <- unicode_clusters <- unicode_features_clusters <- "✘"
  if(length(mean.penalty) > 0){
    #cat("Applying global penalization;\n")
    unicode_global <- "✔"
    Phi <- numeric(R.kept)}
  if("d" %in% mean.penalty){
    #cat("Applying feature penalization;\n")
    unicode_features <- "✔"
    Zeta.d <- matrix(0, R.kept, d)}
  if("c" %in% mean.penalty){
    #cat("Applying cluster penalization;\n")
    unicode_clusters <- "✔"
    Zeta.c <- matrix(0, R.kept, K)}
  if("cd" %in% mean.penalty){
    #cat("Applying feature-cluster penalization;\n")
    unicode_features_clusters <- "✔"
    Zeta.cd <- matrix(0, R.kept, K*d)}
  cat("Penalisation factors:\n")
  cat(paste("global ",unicode_global," |  features  ",unicode_features," |  clusters  ",unicode_clusters," |  features-clusters  ",unicode_features_clusters,sep=""))
  cat("\n")


# Default hyperparameters -------------------------------------------------

  if(is.null(mu0)) mu0 <- rep(0, d)
  if(is.null(Sigma0)){
    if(d > 1) Sigma0 <- diag(rep(10, d)) else
      Sigma0 <- 10}
  shape_tau <- tau.hyperpar[1]
  scale_tau <- tau.hyperpar[2]

# Inizialization ----------------------------------------------------------
  if(!is.null(initialization) & class(initialization) == "perla"){
    max.iter.inizialization <- dim(initialization$Z)[3]
    Z.current <- initialization$Z[,,max.iter.inizialization]
    Rho.current <- initialization$Rho[max.iter.inizialization,]
    Mu.current <- initialization$Mu[,,max.iter.inizialization]
    Sigma.current <- initialization$Sigma[,,max.iter.inizialization]
    Tau.current <- initialization$Tau[max.iter.inizialization]
  } else {   # if the input object is not of class perla, but a simple list...
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
    if(is.null(initialization$Tau))
      Tau.current <- tau.value else Tau.current <- initialization$Tau
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


  detOmega <- determinant((diag(rowSums(W))-rho.value*W)/Tau.current, logarithm=T)$mod
  if(prior.rho == "cont") acceptance.rho <- numeric(K-1) else acceptance.rho <- NULL
  if(prior.rho == "metropolis") acceptance.tau <- 0 else acceptance.tau <- NULL


  pb <- progress_bar$new(total = R) # progress bar


# Covariates --------------------------------------------------------------

  if(!is.null(X)){
    if(class(X) != "data.frame") stop("X is not a 'data.frame' object")
    clean_reg <- lm(y ~ .-1, data = X)
    y <- residuals(clean_reg)
  }


# Cycle -------------------------------------------------------------------

  for(r in 2:R){
    pb$tick()

    # --update mu
    if(is.null(mean.penalty))
      Mu.current <- update_means(y = y, Z = Z.current, Sigma = Sigma.current,
                                 mu0 = mu0, Sigma0 = Sigma0) else {
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
    Z.current <- update_Z(y = y,
                          mu = Mu.current,
                          Sigma = Sigma.current,
                          Psi = Psi.current)

    # --update omega and psi
    updated.psi.omega <- update.psi.omega(psi = Psi.current,
                                          omega = Omega.current,
                                          Z = Z.current,
                                          D = diag(rowSums(W)),
                                          W = W, tau = Tau.current,
                                          rho = Rho.current,
                                          apply.mean.correction = F)
    Psi.current <- updated.psi.omega$psi
    Omega.current <- updated.psi.omega$omega
    if(r %in% tokeep) Prob.current <- convert.to.probabilities(Psi.current)

    # --update tau
    if(prior.tau == "conj"){
      Tau.current <- update_tau(Psi = Psi.current, D = diag(rowSums(W)), W = W, Rho = Rho.current, shape = shape_tau, scale = scale_tau)
    }
    if(prior.tau == "metropolis"){
      Tau.metropolis <- update_tau_metropolis(tau_current = Tau.current, Psi = Psi.current, D = diag(rowSums(W)), W = W, Rho = Rho.current, shape = shape_tau, scale = scale_tau, eps = tau.tuning)
      Tau.current <- Tau.metropolis$value
      acceptance.tau <- acceptance.tau+Tau.metropolis$accepted
    }

    # --update rho
    if(prior.rho == "disc")
      Rho.current <- update_rho_discrete(psi = Psi.current,
                                         tau = Tau.current,
                                         W = W,
                                         lambda.spike = p.spike,
                                         lambda.slab = 1-p.spike,
                                         rho.value = rho.value,
                                         detOmega = detOmega)
    if(prior.rho == "cont"){
      rho.updated <- update_rho_continuous(rho.old = Rho.current,
                                           psi = Psi.current,
                                           tau = Tau.current,
                                           W = W,
                                           lambda_spike = p.spike,
                                           lambda_slab = 1-p.spike, sigma.rho=3)
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
      Tau[r-max(burnin)] <- Tau.current
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
  results <- list(Mu = Mu, Z = Z, Sigma = Sigma, Prob = Prob, Tau = Tau, Rho = Rho,
                  shrinkage.parameters = shrinkage.parameters,
                  acceptance.rho = acceptance.rho, acceptance.tau = acceptance.tau, y = y)
  class(results) <- "perla"
  return(results)

}
