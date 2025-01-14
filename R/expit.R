#' Inverse transformation of the logit function
#'
#' This function implements the inverse transformation of the logit function.
#'
#' @param p a vector of probability values.
#'
#' @return the value `exp(p)/(1+exp(p))`.
#'

expit <- function(p){
  if(any(p > 1 | p < 0)) stop("p contains some values that are not probabilities")
  exp(p)/(1+exp(p))
}
