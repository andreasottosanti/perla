% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{west_states_W}
\alias{west_states_W}
\title{Data on mortality by U.S. counties (neighbor matrix)}
\format{
\subsection{\code{west_states_W}}{

The file contains the neighbor matrix for the data \code{west_states_data}.
}
}
\source{
{WONDER Online Database (https://wonder.cdc.gov/mcd-icd10.html).}
}
\usage{
west_states_W
}
\description{
Data on mortality by U.S. counties from 2016 to 2019 for the West Coast.
Data are available in the WONDER Online Database at
https://wonder.cdc.gov/mcd-icd10.html.
}
\examples{
data(west_states_W)
plot(west_states_W, key=NULL, axis.col = NULL, axis.row = NULL)

}
\keyword{datasets}
