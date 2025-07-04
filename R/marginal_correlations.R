#' Marginal correlation of the outcome variables
#'
#' This function derives the posterior distributions of the marginal correlations across outcome variables based on the MCMC runs.
#'
#' @param values an object of class `perla`.
#'
#' @return an object of class `perla` that contains also `Corr`, an array of the same dimension of `Sigma` collecting the marginal correlation matrix at every MCMC iteration.
#' @export
#'
marginal_correlations <- function(values){
  if(class(values) != "perla") stop("values must be an object of class 'perla'")
  R <- dim(values$Sigma)[3]
  Corr <- values$Sigma
  for(i in 1:R){
    W <- diag(1/sqrt(diag(values$Sigma[,,i])))
    Corr[,,i] <- W %*% Corr[,,i] %*% W
  }
  values$Corr <- Corr
  values
}
