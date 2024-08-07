% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information_criteria.R
\name{information.criteria}
\alias{information.criteria}
\title{Information criteria of \code{perla}}
\usage{
information.criteria(values, loglikelihood.values = NULL)
}
\arguments{
\item{values}{Output of function \code{perla}.}

\item{loglikelihood.values}{The log-likelihood values at each iteration of the
MCMC. If NULL (default) the log-likelihood values are recovered from the
\code{perla} object.}
}
\value{
Returns the information criteria of the \code{perla} object.
}
\description{
This function computes the DIC3 value of a \code{perla} model.
}
