#' `perla` information criteria
#'
#' @description
#' This function computes the DIC3 value of a `perla` model.
#'
#' @param values Output of function `perla`.
#'
#' @param loglikelihood.values The log-likelihood values at each iteration of the
#' MCMC. If NULL (default) the log-likelihood values are recovered directly from the
#' `perla` object passed as input.
#'
#' @return
#' The information criteria of the estimated `perla` model. The smaller is the value, the better is the model fitting.
#'
#' @export
#'

information.criteria <- function(values, loglikelihood.values = NULL){
  if(is.null(loglikelihood.values))
    loglikelihood.values <- recover.loglikelihood(values, type = "m") else
      warning("log-likelihood values passed in input. Are they marginal log-likelihood values?")
  pen <- 2*sum(log(rowMeans(loglikelihood.values$single.val)))
  DIC3 <- -4*mean(loglikelihood.values$val) + pen
  return(DIC3)
}
