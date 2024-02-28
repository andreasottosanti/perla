#' Z from matrix to vector
#' This function convert the selection matrix `Z` used in the model into a vector of length `n` stating to which each observation belong.
#'
#' @param Z an `n` x `K` matrix, or an array of dimension `n` x `K` x `R` (which is the output of the function `perla`).
#'
#' @return A vector of dimension `n` or a matrix of dimension `n` x `R` containing the clustering labels.
#' @export
#'
#' @examples

convert.Z.matrix_to_vector <- function(Z){
  if(is.matrix(Z)) res <- apply(Z, 1, which.max) else
    if (is.array(Z)){
      res <- matrix(0, dim(Z)[3], dim(Z)[1])
      for(r in 1:dim(Z)[3])
        res[r,] <- apply(Z[,,r], 1, which.max)
    }
  res
}

#' Z from vector to matrix
#'
#' This function convert a vector of length `n` stating to which each observation belong into the selection matrix `Z` used in the model
#'
#' @param z
#'
#' @return
#' @export
#'
#' @examples

convert.Z.vector_to_matrix <- function(z){
  Z <- matrix(0, length(z), length(table(z)))
  for(i in 1:nrow(Z)) Z[i, z[i]] <- 1
  Z
}
