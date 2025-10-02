#' Merge multiple `perla` objects into a unique chain
#'
#' @description
#' This function is designed to combine multiple Markov chains obtained running
#' the function `perla` from different starting points.
#'
#' @import abind
#'
#' @param values A list of `perla` objects. These objects must be different
#' parallel runs of the same model.
#'
#' @return An object of class `perla`. It contains also an object called `chains`,
#' which states from which chain the elements come from.
#'
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
    results$Responsabilities <- abind(results$Responsabilities, values[[i]]$Responsabilities)
    results$Sigma <- abind(results$Sigma, values[[i]]$Sigma)
    results$Prob <- abind(results$Prob, values[[i]]$Prob)
    results$Rho <- rbind(results$Rho)
    if("Phi" %in% names(results$shrinkage.parameters))
      results$shrinkage.parameters$Phi <- rbind(results$shrinkage.parameters$Phi, values[[i]]$shrinkage.parameters$Phi)
    if("Zeta.d" %in% names(results$shrinkage.parameters))
      results$shrinkage.parameters$Zeta.d <- rbind(results$shrinkage.parameters$Zeta.d, values[[i]]$shrinkage.parameters$Zeta.d)
    if("Zeta.c" %in% names(results$shrinkage.parameters))
      results$shrinkage.parameters$Zeta.c <- rbind(results$shrinkage.parameters$Zeta.c, values[[i]]$shrinkage.parameters$Zeta.c)
    if("Zeta.cd" %in% names(results$shrinkage.parameters))
      results$shrinkage.parameters$Zeta.cd <- rbind(values[[i]]$shrinkage.parameters$Zeta.cd, results$shrinkage.parameters$Zeta.cd)
    results$acceptance.rho <- rbind(results$acceptance.rho, values[[i]]$acceptance.rho)
    chains <- c(chains, rep(i, dim(values[[i]]$Mu)[3]))
  }
  results$chains <- chains
  class(results) <- "perla"
  return(results)
}
