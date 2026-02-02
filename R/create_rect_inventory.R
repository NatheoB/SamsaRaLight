#' Make the inventory polygon rectangle and optionally axis-align it
#'
#' The minimum-area enclosing rectangle is computed, and if
#' \code{rotate_axisaligned = TRUE}, align the rectangle with the X and Y axes, and 
#' consequently, coordinates of trees (and sensors) are rotated.
#'
#' Tree and sensor coordinates are rotated consistently.
#'
#' @param core_polygon_df A validated data.frame defining the inventory polygon
#'   with columns \code{x} and \code{y}.
#' @param trees A data.frame containing tree coordinates (\code{x}, \code{y}).
#' @param north2x Numeric. Clockwise angle from North to the X-axis (degrees).
#' @param sensors Optional data.frame containing sensor coordinates (\code{x}, \code{y}), or \code{NULL}.
#' @param rotate_axisaligned Logical. If \code{TRUE}, rotate coordinates to align rectangle with axes.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{core_polygon_df}}{Rectangular polygon (rotated axis-aligned if requested).}
#'   \item{\code{trees}}{Updated tree data.frame (rotated if requested).}
#'   \item{\code{sensors}}{Updated sensor data.frame (rotated if requested, or NULL).}
#'   \item{\code{north2x}}{Updated orientation angle (degrees).}
#'   \item{\code{rotation}}{Applied rotation angle in degrees. 0 if no rotation.}
#' }
#'
#' @importFrom sfheaders sf_polygon
#' @importFrom sf st_minimum_rotated_rectangle st_coordinates
#' @importFrom dplyr select distinct
#'
#' @keywords internal
create_rect_inventory <- function(core_polygon_df,
                                  trees,
                                  north2x,
                                  sensors = NULL,
                                  rotate_axisaligned = TRUE) {
  
  ## Convert polygon to sf ----
  core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  
  ## Compute minimum enclosing rectangle ----
  rect_sf <- sf::st_minimum_rotated_rectangle(core_polygon_sf)
  rect_df <- sf::st_coordinates(rect_sf) %>%
    as.data.frame() %>%
    dplyr::select(x = X, y = Y) %>%
    dplyr::distinct()
  
  theta_deg_ccw <- 0  # default if no rotation applied
  
  if (rotate_axisaligned) {
    
    # Determine rotation angle from longest edge with X-axis
    edges <- diff(as.matrix(rbind(rect_df, rect_df[1,]))) 
    lens <- sqrt(rowSums(edges^2)) 
    i <- which.max(lens)
    
    theta_edge <- atan2(edges[[i, "y"]], edges[[i, "x"]])
    
    # nearest multiple of 90 deg
    theta_target <- round(theta_edge / (pi/2)) * (pi/2)
    
    # rotation to apply to align rectangle
    theta <- theta_target - theta_edge
    
    ## Compute centroid
    centroid <- colMeans(rect_df)
    
    ## Rotate rectangle coordinates around the polygon centroid ----
    rect_aligned <- rotate_vec_ccw(rect_df$x - centroid[1], rect_df$y - centroid[2], theta)
    rect_df$x <- rect_aligned$x + centroid[1]
    rect_df$y <- rect_aligned$y + centroid[2]
    
    ## Rotate trees ----
    trees_rot <- rotate_vec_ccw(trees$x - centroid[1], trees$y - centroid[2], theta)
    trees$x <- trees_rot$x + centroid[1]
    trees$y <- trees_rot$y + centroid[2]
    
    ## Rotate sensors if present ----
    if (!is.null(sensors)) {
      sensors_rot <- rotate_vec_ccw(sensors$x - centroid[1], sensors$y - centroid[2], theta)
      sensors$x <- sensors_rot$x + centroid[1]
      sensors$y <- sensors_rot$y + centroid[2]
    }
    
    ## Update north2x (clockwise from north to x axis)
    theta_deg_ccw <- theta * 180 / pi   # math CCW degrees
    north2x <- (north2x + theta_deg_ccw) %% 360 # Thus, add both value because rotation CCW makes the north further of the x-axis
  }
  
  # Update core polygon
  core_polygon_df <- rect_df
  core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  
  ## Return result ----
  return(list(
    core_polygon_df = core_polygon_df,
    core_polygon_sf = core_polygon_sf,
    trees           = trees,
    sensors         = sensors,
    north2x         = north2x,
    rotation_ccw    = theta_deg_ccw
  ))
}

