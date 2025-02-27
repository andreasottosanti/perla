#' Simulate maps using DAGAR
#'
#' @description
#' Simulate maps using DAGAR.
#'
#'
#' @import spdep
#' @import surveillance
#' @import LaplacesDemon
#'
#' @param spatial.map A `SpatialPolygonsDataFrame` object with the map you want
#' to simulate from.
#' @param K Number of clusters (if `z` is passed, it is ignored).
#' @param Sigma (optional) correlation matrix of the variables.
#' @param rho `rho` parameter of the DAGAR.
#' @param d Number of variables (if `Sigma` is passed, it is ignored).
#' @param range.mu Range of the uniform distribution used in data generation.
#' @param range.Sigma Range of the uniform distribution used to generate the
#' correlation matrix Sigma.
#' @param prob.null.centroid Probability that the `K` centroids related to a certain variable are zero. This is used to simulate scenarios where not all the `d` variables are responsible of the formation of the `K` clusters.
#' @param scale.factor.variance Scaling factor for variance-covariance matrix
#' (default 0.05).
#' @param correct.mean.prob Use the stick-breaking formulation of the
#' multinomial distribution while simulating (default `FALSE`).
#' @param plot.map Plot the map during the simulation (default `FALSE`).
#' @param X (optional) a data.frame containing some possible covariates. If passed, each value of the matrix of regression coefficients will be generated using `rnorm(1, 0, range.beta)`.
#' @param z (optional) a vector of clustering labels. If passed, the clusters will not be generated randomly using the multinomial distribution, and `K` will be set equal to the number of unique values of `z`.
#' @param Mu (optional) a matrix of cluster centroids. It must have `K` rows and `d` columns. If passed, the cluster centroids are no longer generated randomly.
#' @param range.beta (optional) the standard deviation of the Gaussian distribution used to generate the regression coefficients
#' @param Beta (optional) a matrix of regression coefficients. It must have `p` rows and `d` columns. If passed, the regression coefficients are no longer generated randomly.
#'
#' @return
#' The function returns the map and the simulated data.
#'
#' @export
#'
#' @examples
#' data(west_states_data)
#'
#' sim1 <- generate.simulations(spatial.map = west_states_data,
#' K = 3,
#' d = 10,
#' Sigma = NULL,
#' range.mu = 0.5,
#' range.Sigma = 0.1,
#' rho = c(0.8, 0.9),
#' prob.null.centroid = 0.5,
#' scale.factor.variance = 0.05,
#' correct.mean.prob = F,
#' plot.map = T)
#'
#'
#' sim2 <- generate.simulations(spatial.map = west_states_data,
#' K = 4,
#' d = 3,
#' Sigma = NULL,
#' range.mu = 1,
#' rho = c(0.01, 0.455, 0.9),
#' range.Sigma = 0.2,
#' prob.null.centroid = 0.5,
#' scale.factor.variance = 0.07,
#' correct.mean.prob = F,
#' plot.map = T)
#'
generate.simulations <- function(spatial.map,
                                 K = NULL,
                                 d = NULL,
                                 X = NULL,
                                 Beta = NULL,
                                 Sigma = NULL,
                                 z = NULL,
                                 Mu = NULL,
                                 range.mu,
                                 range.Sigma,
                                 range.beta = 1,
                                 rho,
                                 prob.null.centroid = 0.5,
                                 scale.factor.variance = 0.05,
                                 correct.mean.prob = F,
                                 plot.map = F
){

  if(!is.null(Sigma)) d <- nrow(Sigma)
  if(!is.null(Mu)) d <- ncol(Mu)
  if(!is.null(Sigma) & !is.null(Mu)){
    if(nrow(Sigma) != ncol(Mu)) stop("The number of rows of Sigma does not match the number of columns of Mu")
  }
  if(!is.null(z)){
    if(any(sort(unique(z)) != (1:length(unique(z))))) stop(paste("z is not a vector of values from 1 to",length(unique(z))))
    K <- max(z)
  }

  # order the map
  new.ordering <- order.map(map = spatial.map, reorder = F)
  W <- new.ordering$W.o
  D <- diag(rowSums(W))

  # define the number of areas
  n <- nrow(new.ordering$map)

  # DAGAR generation: define the probabilities
  map.prob <- NULL
  if(is.null(z)){
    Psi <- matrix(0, n, K-1)
    for(r in 1:(K-1))
      Psi[,r] <- generate.from.DAGAR(rho = rho[r],W.lower = new.ordering$W.o.lower)+rep(-log(K-r), n) * correct.mean.prob
    Prob <- data.frame(convert.to.probabilities(psi = Psi))
    row.names(Prob) <- row.names(new.ordering$map)
    map.prob <- SpatialPolygonsDataFrame(SpatialPolygons(new.ordering$map@polygons), Prob)
  }

  # Deterministic allocation: definition of the clusters using the probabilities
  if(is.null(z)){
    z <- numeric(n)
    z <- apply(Prob, 1, which.max)
  }
  Z <- convert.Z.vector_to_matrix(z)
  if(plot.map) plot(new.ordering$map, col = z)

  # generate Mu
  if(is.null(Mu)){
    Mu <- matrix(0, K, d)
    for(j in 1:d)
      if(runif(1) < (1-prob.null.centroid)) Mu[,j] <- runif(K, -range.mu, range.mu) else Mu[,j] <- 0}

  # generate the regression coefficients, if they aren't given in input
  if(is.null(Beta)){
    if(!is.null(X)){
      if(!is.data.frame(X)) stop("X is not a data.frame object")
      x <- model.matrix(~.-1, X)
      p <- ncol(x)
      B <- matrix(rnorm(p * d, 0, range.beta), p, d)
    }
  } else {
    if(ncol(X) != nrow(Beta)) stop("The number of columns of X does not match the number of rows of Beta")
    B <- Beta
  }


  # generate Sigma, correlation matrix across causes of death, if it isn't given in input
  if(is.null(Sigma)){
    Sigma <- matrix(0, d, d)
    Sigma[lower.tri(Sigma, diag = F)] <- runif(length(Sigma[lower.tri(Sigma, diag = F)]),-range.Sigma,range.Sigma)
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1
    Corr <- round(diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma))),7)
    Sigma <- Corr * scale.factor.variance
  }

  # generate the data
  y <- matrix(0, n, ncol(Sigma))
  eta <- Z %*% Mu
  if(!is.null(X)) eta <- eta + x %*% B
  for(i in 1:n)
    y[i,] <- rmvn(n = 1, mu = eta[i,], Sigma = Sigma)
  row.names(y) <- row.names(new.ordering$map)
  rho.true <- rho
  Mu.true <- Mu
  map <- SpatialPolygonsDataFrame(SpatialPolygons(new.ordering$map@polygons),
                                  data.frame(y))

  # output
  return(list(map = map,
              Mu.true = Mu.true,
              y = y,
              map.prob = map.prob,
              z = z,
              Sigma = Sigma,
              B = B,
              D = D,
              W = W,
              rho.true = rho.true))
}

# function to reorder the spatial map
order.map <- function(map = NULL, coordinates = NULL, W = NULL, reorder = T){
  if(!is.null(map)){
    coordinates <- t(sapply(1:nrow(map), function(i) map@polygons[[i]]@labpt))
    if(is.null(W)) W <- surveillance::poly2adjmat(map)
  }
  if(reorder) o <- order(coordinates[,1], decreasing = F) else o <- 1:nrow(map)
  V <- W[o,o]
  W.lower <- V
  W.lower[upper.tri(W.lower, diag = T)] <- 0
  return(list(order = o, coordinates.o = coordinates[o,], W.o = V,
              W.o.lower = W.lower, map = map[o,]))
}

# function to generate data from a DAGAR model
generate.from.DAGAR <- function(rho, W.lower){
  n <- rowSums(W.lower)
  Tau <- (1+(n-1)*rho^2)/(1-rho^2)
  weights <- rho/(1+(n-1)*rho^2)
  B <- -W.lower * weights
  diag(B) <- 1
  as.vector(LaplacesDemon::rmvnp(1, mu = rep(0, nrow(B)),
                                 Omega = t(B) %*% diag(Tau) %*% B))
}
