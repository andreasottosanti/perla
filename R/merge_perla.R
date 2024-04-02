#' Merge multiple perla objects into a unique chain
#'
#' This function is designed to combine multiple Markov chains obtained running Perla from different starting points.
#' @import abind
#' @param values a list of `perla` objects. They must be different parallel runs of the same model.
#'
#' @return An object of class perla. It contains an object called `chains`, which states from which chain the elements come from.
#' @export
#'
#' @examples

merge.perla <- function(values){
  # check whether all the elements in 'values' are perla objects
  if(any(lapply(values, function(l) class(l)) != "perla")) stop("Not all the input objects are of class perla")

  results <- values[[1]]
  chains <- rep(1, dim(results$Mu)[3])
  for(i in 2:length(values)){
    results$Mu <- abind(results$Mu, values[[i]]$Mu)
    results$Z <- abind(results$Z, values[[i]]$Z)
    results$Sigma <- abind(results$Sigma, values[[i]]$Sigma)
    results$Prob <- abind(results$Prob, values[[i]]$Prob)
    results$Rho <- rbind(results$Rho)
    if("Phi" %in% names(results$shrinkage.parameters)) results$shrinkage.parameters$Phi <- values[[i]]$Phi
    if("Zeta.d" %in% names(results$shrinkage.parameters)) results$shrinkage.parameters$Zeta.d <- values[[i]]$Zeta.d
    if("Zeta.c" %in% names(results$shrinkage.parameters)) results$shrinkage.parameters$Zeta.c <- values[[i]]$Zeta.c
    if("Zeta.cd" %in% names(results$shrinkage.parameters)) results$shrinkage.parameters$Zeta.cd <- values[[i]]$Zeta.cd
    results$acceptance.rho <- rbind(results$acceptance.rho, values[[i]]$acceptance.rho)
    chains <- c(chains, rep(i, dim(values[[i]]$Mu)[3]))
  }
  results$chains <- chains
  class(results) <- "perla"
  return(results)
}
