#' Label switching correction 2
#'
#' This function can be applied to remove the label switching effect from the posterior sample.
#'
#' @import abind
#' @param values an object of class `perla`
#' @param burnin a vector of indexes denoting the MCMC draws to be discarded (default `NULL` means every draw is kept)
#' @param loglikelihood.values the output of the function `loglikelihood.values`. If `NULL` and the methods `PRA` and/or `ECR` are evoced to solve the label switching phenomenon, then the function will automatically compute the loglikelihood values.
#' @param methods default `("PRA","ECR")`
#' @param zpivot if furnished, it performs the relabelling using this vector as reference, otherwise the function will consider the partition corresponding to the largest log-likelihood value
#'
#' @return A list of two elements.
#' - `results` contains the permuted clustering centroids and clustering probabilities for each method used to resolve the label switching;
#' - `relabelling` is the output of the `label.switching` function.
#' @export
#'

remove.label.switching2 <- function(values,
                                   methods = c("PRA","ECR"),
                                   loglikelihood.values = NULL,
                                   burnin = NULL,
                                   zpivot = NULL){
  if(is.null(burnin)) to.keep <- 1:dim(values$Z)[3] else
    to.keep <- setdiff(1:dim(values$Z)[3], burnin)
  if(is.null(loglikelihood.values)){
    #cat("[WARNING]: Log-likelihood values are not provided but required.")
    values <- recover.loglikelihood(values, burnin)
  }
  y <- values$y
  Z <- values$Z[,,to.keep]  # n x K x R
  Mu <- values$Mu[,,to.keep]    # K x d x R
  Prob <- values$Prob[,,to.keep] # n x K x R
  Responsabilities <- values$Responsabilities[,,to.keep]  # n x K x R

  # objects should be transformed as follows:
  # Mu: R x K x d
  # C.values: R x K x 1
  # CD.values: R x K x d
  # Prob: R x K x n
  # and then merged w.r.t. the third dimension

  mcmc.values <- aperm(results$Mu, perm = c(3,1,2))
  if(F){
  if("c" %in% values$mean.penalty){
    C.values <- array(0, dim = c(dim(Z)[3], dim(Z)[2], 1))
    C.values[,,1] <- values$shrinkage.parameters$Zeta.c[to.keep,]  # R x K x 1
    mcmc.values <- abind(mcmc.values, C.values)
  }
  if("cd" %in% values$mean.penalty){
    CD.values <- reshape_Kd_to_array(mat = values$shrinkage.parameters$Zeta.cd[to.keep,], K = dim(Mu)[1], d = dim(Mu)[2])  # R x K x d
    mcmc.values <- abind(mcmc.values, CD.values)
  }
  }
  #mcmc.values <- abind(mcmc.values, aperm(Prob, perm = c(3,2,1)))

  # preparation of the objects ----------------------------------------------
  P <- aperm(Responsabilities, perm = c(3,1,2))   # R x n x K
  Zvec <- convert.Z.matrix_to_vector(Z)
  mapindex <- which.max(values$loglik)


  # relabelling -------------------------------------------------------------
  if(is.null(zpivot)) zpivot <- Zvec[mapindex,]
  relabelling <- label.switching(method = methods, data = y,
                                 mcmc = mcmc.values, prapivot = mcmc.values[mapindex,,],
                                 z = Zvec, zpivot = zpivot, p = P, K = dim(Mu)[1])

  results <- list()
  for(m in 1:length(methods)){
    Permutations <- permute.mcmc(mcmc = mcmc.values, permutations = relabelling$permutations[[m]])$out
    results[[m]] <- list()
    results[[m]]$M <- Permutations[,,1:(dim(Mu)[2])]
    if(F){
    if("c" %in% values$mean.penalty){
      results[[m]]$C.pen <- Permutations[,,(dim(Mu)[2])+1]
    }
    if("cd" %in% values$mean.penalty){
      if("c" %in% values$mean.penalty){
        results[[m]]$CD.pen <- Permutations[,,(dim(Mu)[2]+2):(2*dim(Mu)[2]+1)]
      } else {
        results[[m]]$CD.pen <- Permutations[,,(dim(Mu)[2]+1):(2*dim(Mu)[2])]
      }
    }
    }
    #results[[m]]$P <- Permutations[,,(dim(Permutations)[3]-nrow(y)+1):dim(Permutations)[3]]
    results[[m]]$Responsabilities <- permute.mcmc(mcmc = aperm(Responsabilities, perm = c(3,2,1)), permutations = relabelling$permutations[[m]])$out
  }
  names(results) <- methods
  return(list(results = results, relabelling = relabelling))
}




reshape_Kd_to_array <- function(mat, K, d) {
  # mat: matrix of dimension n x (K*d)
  # K: number of clusters
  # d: dimension of each centroid

  n <- nrow(mat)
  array_data <- array(NA, dim = c(n, K, d))

  for (i in 1:d) {
    # Extract columns for dimension i (each dimension occupies K columns)
    start_col <- (i - 1) * K + 1
    end_col <- i * K
    array_data[,,i] <- mat[, start_col:end_col, drop = FALSE]
  }

  return(array_data)
}



