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
#' @param K Number of clusters.
#' @param Sigma Correlation matrix across the events that are inducing the
#' clusters.
#' @param rho `rho` parameter of the DAGAR.
#' @param d Number of the events that are inducing the clusters.
#' @param range.mu Range of the uniform distribution used in data generation.
#' @param prob.null.centroid Probability of the uniform distribution used in
#' data generation.
#' @param scale.factor.variance Scaling factor for variance-covariance matrix (
#' default 0.05).
#' @param correct.mean.prob Use the stick-breaking formulation of the
#' multinomial distribution while simulating (default `FALSE`).
#' @param plot.map Plot the map during the simulation (default `FALSE`).
#'
#' @return
#' The function returns the map and the simulated data.
#'
#' @export
#'
#' @examples
#' west.states <- data(west_states_data)
#'
#' sim1 <- generate.simulations(spatial.map = west.states,
#' K = 3,
#' d = 10,
#' Sigma = NULL,
#' range.mu = 0.5,
#' rho = c(0.8, 0.9),
#' prob.null.centroid = 0.5,
#' scale.factor.variance = 0.05,
#' correct.mean.prob = F,
#' plot.map = F)
#'
#'Corr.matrix <- matrix(c(0.07, -0.01175893, 0.00124239,
#'-0.01175893, 0.07, 0.00653324,
#'0.00124239, 0.00653324, 0.07),
#'byrow = T, ncol = 3, nrow = 3)
#'
#' sim2 <- generate.simulations(spatial.map = west.states,
#' K = 3,
#' d = NULL,
#' Sigma = Corr.matrix,
#' range.mu = 1, rho = c(0.01, 0.9),
#' prob.null.centroid = 0.5,
#' scale.factor.variance = 0.05,
#' correct.mean.prob = F,
#' plot.map = F)
#'
generate.simulations <- function(spatial.map,
                                 K,
                                 d = NULL,
                                 Sigma = NULL,
                                 range.mu,
                                 rho,
                                 prob.null.centroid = 0.5,
                                 scale.factor.variance = 0.05,
                                 correct.mean.prob = F,
                                 plot.map = F
){

  if(!is.null(Sigma)) d <- nrow(Sigma)

  # order the map
  new.ordering <- order.map(map = spatial.map, reorder = F)
  W <- new.ordering$W.o
  D <- diag(rowSums(W))

  # define the number of areas
  n <- nrow(new.ordering$map)

  # DAGAR generation: define the probabilities
  Psi <- matrix(0, n, K-1)
  for(r in 1:(K-1))
    Psi[,r] <- generate.from.DAGAR(rho = rho[r],W.lower = new.ordering$W.o.lower)+rep(-log(K-r), n) * correct.mean.prob
  Prob <- data.frame(convert.to.probabilities(psi = Psi))
  row.names(Prob) <- row.names(new.ordering$map)
  map.prob <- SpatialPolygonsDataFrame(SpatialPolygons(new.ordering$map@polygons), Prob)

  # Deterministic allocation: definition of the clusters using the probabilities
  z <- numeric(n)
  z <- apply(Prob, 1, which.max)
  Z <- convert.Z.vector_to_matrix(z)
  if(plot.map) plot(new.ordering$map, col = z)

  # generate Mu
  Mu <- matrix(0, K, d)
  for(j in 1:d)
    if(runif(1) < (1-prob.null.centroid)) Mu[,j] <- runif(K, -range.mu, range.mu) else Mu[,j] <- 0

  # generate Sigma, correlation matrix across causes of death, if it isn't given in input
  if(is.null(Sigma)){
    Sigma <- matrix(0, d, d)
    Sigma[lower.tri(Sigma, diag = F)] <- runif(length(Sigma[lower.tri(Sigma, diag = F)]),-.1,.1)
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1
    Corr <- round(diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma))),7)
    Sigma <- Corr * scale.factor.variance
  }

  # generate the data
  y <- matrix(0, n, ncol(Sigma))
  eta <- Z %*% Mu
  for(i in 1:n)
    y[i,] <- rmvn(n = 1, mu = eta[i,], Sigma = Sigma)
  row.names(y) <- row.names(Prob)
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
