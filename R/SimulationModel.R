#' Simulate maps
#'
#' @import spdep
#' @import surveillance
#'
#' @param map
#' @param K
#' @param Sigma
#' @param rho
#' @param simulation.model
#' @param tau
#' @param allocation
#' @param d
#' @param range.mu
#' @param prob.null.centroid
#' @param scale.factor.variance
#' @param correct.mean.prob
#' @param plot.map
#'
#' @return
#' @export
#'
#' @examples
generate.simulations <- function(map,
                                 K,
                                 d = NULL,
                                 Sigma = NULL,
                                 range.mu,
                                 rho,
                                 simulation.model = "dagar",
                                 tau = NULL,
                                 allocation = "deter",
                                 prob.null.centroid = 0.5,
                                 scale.factor.variance = 0.005,
                                 correct.mean.prob = F,
                                 plot.map = F
){
  if(!is.null(Sigma)) d <- nrow(Sigma)
  new.ordering <- order.map(map = map, reorder = F)
  W <- new.ordering$W.o
  D <- diag(rowSums(W))
  n <- nrow(new.ordering$map)
  z <- numeric(n)
  Psi <- matrix(0, n, K-1)
  if(is.null(tau)) tau <- rep(1, K-1)

  # CAR generation ----------------------------------------------------------
  if(simulation.model == "car"){
    for(r in 1:(K-1))
      Psi[,r] <- rmvnp(n = 1, mu = rep(0,n), Omega = (D - rho[r]*W)/tau[1])+
        rep(-log(K-r), n) * I(correct.mean.prob)
  }

  # DAGAR generation --------------------------------------------------------
  if(simulation.model == "dagar"){
    for(r in 1:(K-1))
      Psi[,r] <- generate.from.DAGAR(rho = rho[r], W.lower = new.ordering$W.o.lower)+
        rep(-log(K-r), n) * correct.mean.prob
  }

  Prob <- data.frame(convert.to.probabilities(psi = Psi))
  row.names(Prob) <- row.names(new.ordering$map)
  map.prob <- SpatialPolygonsDataFrame(SpatialPolygons(new.ordering$map@polygons), Prob)


# Allocation --------------------------------------------------------------
  # Probabilistic allocation ------------------------------------------------
  if(allocation == "prob"){
    Z <- t(apply(Prob, 1, function(p) rmultinom(n = 1, size = 1, prob = p)))
    z <- Z.fromMatrix_toVector(Z)
    if(plot.map) plot(new.ordering$map, col = z)
  }

  # Deterministic allocation ------------------------------------------------
  if(allocation == "deter"){
    z <- apply(Prob, 1, which.max)
    Z <- convert.Z.vector_to_matrix(z)
    if(plot.map) plot(new.ordering$map, col = z)
  }

  Mu <- matrix(0, K, d)
  #for(k in 1:3)
  for(j in 1:d)
    if(runif(1) < (1-prob.null.centroid)) Mu[,j] <- runif(K, -range.mu, range.mu) else Mu[,j] <- 0
  if(is.null(Sigma)){
    Sigma <- matrix(0, d, d)
    Sigma[lower.tri(Sigma, diag = F)] <- runif(length(Sigma[lower.tri(Sigma, diag = F)]),-.1,.1)
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- 1
    Corr <- round(diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma))),7)
    Sigma <- Corr * scale.factor.variance
  }
  y <- matrix(0, n, ncol(Sigma))
  eta <- Z %*% Mu
  for(i in 1:n)
    y[i,] <- rmvn(n = 1, mu = eta[i,], Sigma = Sigma)
  row.names(y) <- row.names(Prob)
  rho.true <- rho
  Mu.true <- Mu
  map <- SpatialPolygonsDataFrame(SpatialPolygons(new.ordering$map@polygons), data.frame(y))
  # Y <- data.frame(y = as.vector(exp(y)),
  #                 cause = (paste("cause",sapply(1:d, function(i) rep(i, nrow(y))))),
  #                 cluster = rep(factor(z),d))
  #Y$cause <- factor(Y$cause, levels = unique(Y$cause))

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


order.map <- function(map = NULL, coordinates = NULL, W = NULL, reorder = T){
  if(!is.null(map)){
    coordinates <- t(sapply(1:nrow(map), function(i) map@polygons[[i]]@labpt))
    if(is.null(W)) W <- poly2adjmat(map)
  }
  if(reorder) o <- order(coordinates[,1], decreasing = F) else o <- 1:nrow(map)
  V <- W[o,o]
  W.lower <- V
  W.lower[upper.tri(W.lower, diag = T)] <- 0
  return(list(order = o, coordinates.o = coordinates[o,], W.o = V, W.o.lower = W.lower, map = map[o,]))
}

generate.from.DAGAR <- function(rho, W.lower){
  n <- rowSums(W.lower)
  Tau <- (1+(n-1)*rho^2)/(1-rho^2)
  weights <- rho/(1+(n-1)*rho^2)
  B <- -W.lower * weights
  diag(B) <- 1
  as.vector(rmvnp(1, mu = rep(0, nrow(B)), Omega = t(B) %*% diag(Tau) %*% B))
}
