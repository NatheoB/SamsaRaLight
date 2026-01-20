#' Axis-align an inventory polygon and rotate coordinates accordingly
#'
#' This internal function ensures that the inventory zone is an axis-aligned
#' rectangle. If the input polygon is not already a rectangle, the minimum-area
#' enclosing rectangle is computed and rotated to align with the X and Y axes
#' using the minimum possible rotation.
#'
#' Tree and sensor coordinates are rotated consistently.
#'
#' @param core_polygon_df A validated data.frame defining the inventory polygon
#'   with columns \code{x} and \code{y}
#' @param trees_inv A data.frame containing tree coordinates (\code{x}, \code{y})
#' @param north2x Numeric. Clockwise angle from North to the X-axis (degrees)
#' @param sensors Optional data.frame containing sensor coordinates
#'   (\code{x}, \code{y}), or \code{NULL}
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{core_polygon_df}}{Axis-aligned rectangular polygon}
#'   \item{\code{trees_inv}}{Updated tree data.frame}
#'   \item{\code{sensors}}{Updated sensor data.frame (or NULL)}
#'   \item{\code{north2x}}{Updated orientation angle}
#'   \item{\code{rotation_rad}}{Applied rotation angle (radians)}
#' }
#'
#' @keywords internal
create_aarect_inventory <- function(core_polygon_df,
                                    trees_inv,
                                    north2x,
                                    sensors = NULL) {
  
  
  ## Helper: rotate coordinates ----
  rotate_xy <- function(x, y, theta, cx, cy) {
    x0 <- x - cx
    y0 <- y - cy
    xr <-  x0 * cos(theta) + y0 * sin(theta)
    yr <- -x0 * sin(theta) + y0 * cos(theta)
    cbind(x = xr + cx, y = yr + cy)
  }
  
  
  ## Convert polygon to sf ----
  core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  
  
  ## Create minimum enclosing rectangle only if needed ----
  rect_sf <- sf::st_minimum_rotated_rectangle(core_polygon_sf)

  
  ## Extract rectangle coordinates ----
  rect_df <- sf::st_coordinates(rect_sf) %>%
    as.data.frame() %>%
    dplyr::select(x = X, y = Y) %>%
    dplyr::distinct()
  
  
  ## Determine rotation angle from longest edge ----
  edges <- diff(as.matrix(rbind(rect_df, rect_df[1, ])))
  lens  <- sqrt(rowSums(edges^2))
  i     <- which.max(lens)
  
  theta_raw <- atan2(edges[i, "y"], edges[i, "x"])[[1]]
  
  
  ## Apply minimum rotation to reach axis alignment ----
  # Find closest multiple of 90° (pi/2)
  k <- round(theta_raw / (pi / 2))
  theta <- theta_raw - k * (pi / 2)
  
  
  ## Rotate rectangle to align with axes ----
  centroid <- colMeans(rect_df)
  
  rect_aligned <- rotate_xy(
    rect_df$x, rect_df$y,
    theta = theta,
    cx = centroid[1],
    cy = centroid[2]
  ) %>% as.data.frame()
  
  core_polygon_df <- rect_aligned
  core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  
  
  ## Rotate trees ----
  trees_xy <- rotate_xy(
    trees_inv$x, trees_inv$y,
    theta = theta,
    cx = centroid[1],
    cy = centroid[2]
  )
  trees_inv$x <- trees_xy[, "x"]
  trees_inv$y <- trees_xy[, "y"]
  
  
  ## Rotate sensors if present ----
  if (!is.null(sensors)) {
    sensors_xy <- rotate_xy(
      sensors$x, sensors$y,
      theta = theta,
      cx = centroid[1],
      cy = centroid[2]
    )
    sensors$x <- sensors_xy[, "x"]
    sensors$y <- sensors_xy[, "y"]
  }
  
  
  ## Update north2x ----
  theta_deg <- rad2deg(theta)
  north2x <- (north2x - theta_deg) %% 360
  
  
  ## Output ----
  return(list(
    core_polygon_df = core_polygon_df,
    trees_inv       = trees_inv,
    sensors         = sensors,
    north2x         = north2x,
    rotation        = theta_deg
  ))
}
