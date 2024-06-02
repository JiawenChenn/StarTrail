#' Function to keep the points that are within to the boarder of a given coordinates
#'
#' @param coord_mesh a matrix of new coordinates (e.g. super-resolution grid)
#' @param coords a matrix of original coordinates 
#' @param close_dis the distance to the boundary for a point to be considered as a boundary point
#' @param distance_to_boundary add some extra distance to the boundary (default is 0)
#' @returns a logical vector indicating whether the point is kept
#'

coord_filter <- function(coord_mesh, coords, close_dis = 0.2, distance_to_boundary = 0) {
  
  # Ensure the inputs are matrices
  coord_mesh <- as.matrix(coord_mesh)
  coords <- as.matrix(coords)
  
  n <- nrow(coord_mesh)
  indices <- logical(n)
  
  for (i in 1:n) {
    close_point_indices_x <- which(abs(coords[, 1] - coord_mesh[i, 1]) < close_dis)
    close_point_indices_y <- which(abs(coords[, 2] - coord_mesh[i, 2]) < close_dis)
    
    if (length(close_point_indices_x) < 2 || length(close_point_indices_y) < 2) {
      indices[i] <- FALSE
    } else {
      close_point_x <- coords[close_point_indices_x, ]
      close_point_y <- coords[close_point_indices_y, ]
      
      in_y_range <- (coord_mesh[i, 2] < max(close_point_x[, 2]) - distance_to_boundary) & 
                    (coord_mesh[i, 2] > min(close_point_x[, 2]) + distance_to_boundary)
      
      in_x_range <- (coord_mesh[i, 1] < max(close_point_y[, 1]) - distance_to_boundary) & 
                    (coord_mesh[i, 1] > min(close_point_y[, 1]) + distance_to_boundary)
      
      indices[i] <- in_y_range & in_x_range
    }
  }
  
  return(indices)
}