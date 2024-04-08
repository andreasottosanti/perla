#' A posteriori log-likelihood computation
#'
#' This function returns the log-likelihood value at each iteration of the MCMC.
#'
#' @param values an object of class `perla`
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#' @param type if `j`, it returns the joint `(y,z)` log-likelihood values; if `m`, it returns the marginal log-likelihood values of `y`.
#'
#' @return
#' @export
#'
#' @examples
#'

recover.loglikelihood <- function(values, burnin = NULL, type = "j"){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  y <- values$y
  z <- values$Z[,,to.keep]  # n x K x R
  Mu <- values$Mu[,,to.keep]    # K x d x R
  Prob <- values$Prob[,,to.keep]  # n x K x R
  Sigma <- values$Sigma[,,to.keep] # d x d x R
  n <- nrow(y)
  d <- ncol(y)
  K <- dim(Mu)[1]
  val <- numeric(dim(values$Z)[3])
  pb <- progress_bar$new(total = dim(values$Z)[3])
  if(type == "j"){
    for(r in 1:dim(values$Z)[3]){
      pb$tick()
      eta <- values$Z[,,r] %*% Mu[,,r]
      val[r] <- -.5*sum(diag((y - eta) %*% solve(Sigma[,,r]) %*% t(y - eta)))-
        n*d/2*log(2*pi)-n/2*determinant(Sigma[,,r], logarithm = T)$mod
    }
    output <- val
  } else if(type == "m"){
    single.values <- matrix(0, n, dim(values$Z)[3])
      for(r in 1:dim(values$Z)[3]){
        pb$tick()
        likelihood.values.r <- numeric(nrow(y))
        for(k in 1:K){
          eta <- matrix(Mu[k,,r], n, d, byrow = T)
          likelihood.values.r <- likelihood.values.r +
            exp(-.5*(diag((y - eta) %*% solve(Sigma[,,r]) %*% t(y - eta)))-
                  d/2*log(2*pi)-.5*determinant(Sigma[,,r], logarithm = T)$mod+
                  log(Prob[,k,r]))
        }
        single.values[,r] <- likelihood.values.r
        val[r] <- sum(log(likelihood.values.r))
      }
    output <- list(val = val, single.val = single.values)
    } else
      stop("log-likelihood type not recognized. Please select 'm' or 'j'")
  output
}
