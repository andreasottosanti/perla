#' Data on mortality by U.S. counties (`SpatialPolygonDataFrame` object)
#'
#' Data on mortality by U.S. counties from 2016 to 2019 for the West Coast.
#' Data are available in the WONDER Online Database at
#' https://wonder.cdc.gov/mcd-icd10.html.
#'
#'
#' @format ## `west_states_data`
#' The file contains a `SpatialPolygonDataFrame` object with the map and the
#' data of the mortality by U.S. counties from 2016 to 2019 for three causes
#' of death (circulatory diseases, respiratory diseases and tumors).
#'
#' @source {WONDER Online Database (https://wonder.cdc.gov/mcd-icd10.html).}
#'
#' @import sp
#'
#' @examples
#' library(sp)
#' data(west_states_data)
#' plot(west_states_data)
#'
"west_states_data"


#' Data on mortality by U.S. counties (neighbor matrix)
#'
#' Data on mortality by U.S. counties from 2016 to 2019 for the West Coast.
#' Data are available in the WONDER Online Database at
#' https://wonder.cdc.gov/mcd-icd10.html.
#'
#'
#' @format ## `west_states_W`
#' The file contains the neighbor matrix for the data `west_states_data`.
#'
#' @source {WONDER Online Database (https://wonder.cdc.gov/mcd-icd10.html).}
#'
#' @import sp
#' @import plot.matrix
#'
#' @examples
#' data(west_states_W)
#' plot(west_states_W, key=NULL, axis.col = NULL, axis.row = NULL)
#'
"west_states_W"
