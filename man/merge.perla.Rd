% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_perla.R
\name{merge.perla}
\alias{merge.perla}
\title{Merge multiple \code{perla} objects into a unique chain}
\usage{
merge.perla(values)
}
\arguments{
\item{values}{A list of \code{perla} objects. These objects must be different
parallel runs of the same model.}
}
\value{
An object of class \code{perla}. It contains an object called \code{chains},
which states from which chain the elements come from.
}
\description{
This function is designed to combine multiple Markov chains obtained running
the function \code{perla} from different starting points.
}
