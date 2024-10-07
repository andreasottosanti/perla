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
#' Returns the information criteria of the estimated `perla` model.
#'
#' @export
#'
#' @examples
#'

information.criteria <- function(values, loglikelihood.values = NULL){
  if(is.null(loglikelihood.values))
    loglikelihood.values <- recover.loglikelihood(values, type = "m") else
      warning("log-likelihood values passed in input. Are they marginal log-likelihood values?")
  ll <- -2*loglikelihood.values$val
  pen <- 2*sum(log(colMeans(loglikelihood.values$single.val)))
  DIC3 <- 2*mean(ll) + pen
  return(DIC3)
}
