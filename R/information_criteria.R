#' Title
#'
#' @param values an object of class `perla`
#' @param loglikelihood.values
#'
#' @return
#' @export
#'
#' @examples
information.criteria <- function(values, loglikelihood.values = NULL){
  if(is.null(loglikelihood.values))
    loglikelihood.values <- recover.loglikelihood(values)
  ll <- -2*loglikelihood.values
  Ez <- matrix(0, nrow(values))
  DIC4 <- 2*mean(ll)

}
