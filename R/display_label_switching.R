#' Display label switching correction
#'
#' This function displays the relabelling of the cluster centroids $\mu_{k}$ based on the applied algorithms.
#' @import ggplot2
#' @param values the output of the function `remove.label.switching`
#' @param colors (optional) colors of the clusters in the traceplots
#' @param names.diseases (optional) the disease names
#'
#' @returns a `ggplot` object
#' @export
#'
#' @examples
display.label.switching <- function(values, colors = NULL, names.diseases = NULL){
  reshape <- reshape_relabelling(values$results)
  if(!is.null(names.diseases)) levels(reshape$disease) <- names.diseases
  p <- ggplot(reshape, aes(x = iteration, y = value, col = cluster))+
    geom_line(alpha = 0.7)+
    facet_grid(method ~ disease,
               scales="free_y",
               switch = "y"
    )+
    theme_bw()+
    labs(x = "Iteration", y = expression(mu), col = "Cluster")+
    theme(text = element_text(size = 16),
          legend.position = "bottom")
    if(!is.null(colors)) p <- p + scale_color_manual(values = colors)
  p
}



reshape_relabelling <- function(x) {
  result_list <- lapply(names(x), function(method_name) {
    M <- x[[method_name]]$M
    R <- dim(M)[1]
    K <- dim(M)[2]
    d <- dim(M)[3]

    # Create a data frame by expanding the grid
    df <- do.call(rbind, lapply(1:K, function(k) {
      do.call(rbind, lapply(1:d, function(j) {
        data.frame(
          iteration = 1:R,
          value = M[, k, j],
          cluster = k,
          disease = j,
          method = method_name
        )
      }))
    }))

    return(df)
  })

  # Combine all method-specific data frames into one
  final_df <- do.call(rbind, result_list)
  rownames(final_df) <- NULL
  final_df$cluster <- factor(final_df$cluster)
  final_df$disease <- factor(final_df$disease)
  final_df$method <- factor(final_df$method)
  return(final_df)
}
