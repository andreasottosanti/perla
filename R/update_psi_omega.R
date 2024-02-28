update.psi.ik <- function(psi_k,
                          Z_ik,
                          N_ik,
                          omega_ik,
                          D_ii,
                          W_i,
                          tau,
                          rho = .99){
  v_ik <- Z_ik - N_ik/2
  u_ik <- rho * sum(W_i * psi_k)/D_ii
  mean.comp <- v_ik + u_ik * D_ii / tau
  var.comp <- D_ii / tau + omega_ik
  rnorm(1, mean = mean.comp/var.comp, sqrt(1/var.comp))
}


update.psi.omega <- function(psi, omega, Z, D, W, tau, rho = 0.99, apply.mean.correction = F){
  n <- nrow(Z)
  if(length(rho) == 1 & !is.vector(psi)) rho <- rep(rho, ncol(psi))
  if(is.vector(psi)) psi <- matrix(psi, ncol = 1)
  if(is.vector(omega)) omega <- matrix(omega, ncol = 1)
  for(i in 1:n){
    for(k in 1:ncol(psi)){
      N_ik <- ifelse(k > 1, 1 - sum(Z[i,1:(k-1)]), 1)
      omega[i,k] <- ifelse(N_ik > 0, pgdraw(N_ik, psi[i,k]), 0)
      psi[i,k] <- -log(ncol(psi)+1-k)*apply.mean.correction+
        update.psi.ik(psi_k = psi[,k], Z_ik = Z[i,k],
                                N_ik = N_ik, omega_ik = omega[i,k],
                                D_ii = D[i,i], W_i = W[i,],
                                tau = tau, rho = rho[k])
    }
  }
  list(psi = psi, omega = omega)
}
