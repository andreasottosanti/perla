#' `perla` information criteria
#'
#' @description
#' This function computes the DIC3 value of a `perla` model.
#'
#' @param values an object of class `perla`.
#'
#' @param loglikelihood.values The log-likelihood values at each iteration of the
#' MCMC. If NULL (default) the log-likelihood values are recovered directly from the
#' `perla` object passed as input.
#'
#' @return
#' an object of class `perla` that contains also `DIC3`, the information criterion value of the estimated `perla` model. The smaller is the value, the better is the model fitting.
#'
#' @export
#'

information.criteria <- function(values, loglikelihood.values = NULL){
  if(class(values) != "perla") stop("values must be an object of class 'perla'")
  if(is.null(loglikelihood.values))
    loglikelihood.values <- recover.loglikelihood.marginal(values) else
      warning("log-likelihood values passed in input. Are they marginal log-likelihood values?")
  pen <- 2*sum(log(rowMeans(loglikelihood.values$single.val)))
  values$DIC3 <- -4*mean(loglikelihood.values$val) + pen
  return(values)
}
