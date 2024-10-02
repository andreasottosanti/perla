
update_rho_discrete <- function(psi,
                             tau = 1, W, lambda.spike = 1, lambda.slab = 1, rho.value = .99, detOmega = NULL){
  D <- rowSums(W)
  if(is.vector(psi)) psi <- matrix(psi, ncol = 1)
  rho <- numeric(ncol(psi))
  Omega <- (diag(D)-rho.value*W)/tau
  if(is.null(detOmega)) detOmega <- determinant(Omega, logarithm = T)$mod
  kernelGaussian <- t(psi) %*% Omega %*% psi
  for(k in 1:ncol(psi)){
    log.prob.spike <- sum(dnorm(psi[,k], 0, sqrt(tau/D), log = T))
    log.prob.slab <- ((detOmega) - kernelGaussian[k,k])/2-length(D)/2*log(2*pi)
    prob.spike <- lambda.spike/(lambda.spike+lambda.slab) * exp(log.prob.spike)
    prob.slab <- lambda.slab/(lambda.spike+lambda.slab) * exp(log.prob.slab)
    rho[k] <- sample(0:1, size = 1, prob = c(prob.spike, prob.slab))
    if(rho[k] == 1) rho[k] <- rho.value}
  return(rho)
}
