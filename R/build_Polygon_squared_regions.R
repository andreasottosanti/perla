#' Build a `SpatialPolygon` object made of squared regions
#'
#' This function creates a rectangular `SpatialPolygons` map divided into square regions. The user specifies the number of rows and columns to define the grid's layout.
#'
#' @import sp
#' @param n_rows the number of rows of the map
#' @param n_cols the number of columns of the map
#'
#' @return A `SpatialPolygon` object.
#' @export
#'
#' @examples
#' build_Polygon_squared_regions(n_rows = 5, n_cols = 4)
#'
build_Polygon_squared_regions <- function(n_rows = 5, # Number of rows
                                          n_cols = 4 # Number of columns
                                          )
  {


  square_size <- 1 # Size of each square (length of one side)

  # Create a list to store the polygons
  polygons <- list()

  # Loop to create squares and store them as polygons
  id <- 1
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      x_min <- (j - 1) * square_size
      x_max <- j * square_size
      y_min <- (n_rows - i) * square_size
      y_max <- (n_rows - i + 1) * square_size

      # Define the corners of the square
      coords <- matrix(c(
        x_min, y_min,
        x_min, y_max,
        x_max, y_max,
        x_max, y_min,
        x_min, y_min
      ), ncol = 2, byrow = TRUE)

      # Create a Polygon object
      polygon <- Polygon(coords)
      polygons[[id]] <- Polygons(list(polygon), as.character(id))
      id <- id + 1
    }
  }

  # Create a SpatialPolygons object
  spatial_polygons <- SpatialPolygons(polygons)

  spatial_polygons
}
