#' Label switching correction
#'
#' This function can be applied to remove the label switching effect from the posterior sample.
#'
#' @param values an object of class `perla`
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#' @param loglikelihood.values the output of the function `loglikelihood.values`. If `NULL` and the methods `PRA` and/or `ECR` are evoced to solve the label switching phenomenon, then the function will automatically compute the loglikelihood values.
#' @param methods default `("PRA","ECR")`
#' @param zpivot if furnished, it performs the relabelling using this vector as reference, otherwise the function will consider the partition corresponding to the largest log-likelihood value
#'
#' @return a list of two elements.
#' - `results` contains the permuted clustering centroids and clustering probabilities for each method used to resolve the label switching;
#' - `relabelling` is the output of the `label.switching` function.
#' @export
#'

remove.label.switching <- function(values,
                                  methods = c("PRA","ECR"),
                                  loglikelihood.values = NULL, burnin = NULL, zpivot = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  if(is.null(loglikelihood.values)){
    #cat("[WARNING]: Log-likelihood values are not provided but required.")
    loglikelihood.values <- recover.loglikelihood(values, burnin)
  }
  y <- values$y
  Z <- values$Z[,,to.keep]  # n x K x R
  Mu <- values$Mu[,,to.keep]    # K x d x R
  Prob <- values$Prob[,,to.keep] # n x K x R

# preparation of the objects ----------------------------------------------
  M <- aperm(Mu, perm = c(3,1,2))     # R x K x d
  P <- aperm(Prob, perm = c(3,1,2))   # R x n x K
  Zvec <- convert.Z.matrix_to_vector(Z)
  mapindex <- which.max(loglikelihood.values)


# relabelling -------------------------------------------------------------
  if(is.null(zpivot)) zpivot <- Zvec[mapindex,]
  relabelling <- label.switching(method = methods,
                                 mcmc = M, prapivot = M[mapindex,,],
                                 z = Zvec, zpivot = zpivot, p = P, K = dim(Mu)[1])

  results <- list()
  for(m in 1:length(methods)){
    results[[m]] <- list()
    results[[m]]$M <- permute.mcmc(mcmc = M, permutations = relabelling$permutations[[m]])$out
    results[[m]]$P <- permute.mcmc(mcmc = aperm(Prob, perm = c(3,2,1)), permutations = relabelling$permutations[[m]])$out
  }
  names(results) <- methods
  return(list(results = results, relabelling = relabelling))
}
