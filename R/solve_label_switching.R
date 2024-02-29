#' Solve label switching
#'
#' This function can be applied to remove the label switching effect.
#'
#' @param values an object of class `perla`
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#' @param loglikelihood.values
#' @param methods
#'
#' @return
#' @export
#'
#' @examples

solve.label.switching <- function(values,
                                  methods = c("PRA","ECR","STEPHENS"),
                                  loglikelihood.values = NULL, burnin = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  if(is.null(loglikelihood.values)) loglikelihood.values <- model.loglikelihood(values, burnin)
  y <- values$y
  Z <- values$Z[,,to.keep]  # n x K x R
  Mu <- values$Mu[,,to.keep]    # K x d x R
  Prob <- values$Prob[,,to.keep] # n x K x R

# preparation of the objects ----------------------------------------------
  M <- aperm(Mu, perm = c(3,1,2))
  P <- aperm(Prob, perm = c(3,1,2))
  Zvec <- convert.Z.matrix_to_vector(Z)
  mapindex <- which.max(loglikelihood.values)


# relabelling -------------------------------------------------------------
  relabelling <- label.switching(method = methods,
                                 mcmc = M, prapivot = M[mapindex,,],
                                 z = Zvec, zpivot = Zvec[mapindex,], p = P)

  results <- list()
  for(m in 1:length(methods)){
    results[[m]] <- list()
    results[[m]]$M <- permute.mcmc(mcmc = M, permutations = relabelling$permutations$PRA)

  }
  names(results) <- methods
  relabelling$similarity
  Mu.perm.PRA <-
  Mu.perm.ECR <- permute.mcmc(mcmc = M, permutations = relabelling$permutations$ECR)
  mean.clusters.PRA <- mcmc.list(mcmc((Mu.perm.PRA$output[,1,])),  # elements of cluster 1
                                 mcmc((Mu.perm.PRA$output[,2,])),
                                 mcmc((Mu.perm.PRA$output[,3,])))#,
  #mcmc((Mu.perm.PRA$output[,4,])))
  mean.clusters.ECR <- mcmc.list(mcmc((Mu.perm.ECR$output[,1,])),  # elements of cluster 1
                                 mcmc((Mu.perm.ECR$output[,2,])),
                                 mcmc((Mu.perm.ECR$output[,3,])))#,
  #mcmc((Mu.perm.ECR$output[,4,])))

}
