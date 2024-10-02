convert.Z.matrix_to_vector <- function(Z){
  if(is.matrix(Z)) res <- apply(Z, 1, which.max) else
    if (is.array(Z)){
      res <- matrix(0, dim(Z)[3], dim(Z)[1])
      for(r in 1:dim(Z)[3])
        res[r,] <- apply(Z[,,r], 1, which.max)
    }
  res
}


convert.Z.vector_to_matrix <- function(z){
  Z <- matrix(0, length(z), length(table(z)))
  for(i in 1:nrow(Z)) Z[i, z[i]] <- 1
  Z
}
