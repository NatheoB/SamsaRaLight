# FILE FOR ALL UTILS FUNCTIONS #

#' Convert radians into degrees
#' @noRd
rad2deg <- function(rad) {(rad * 180) / pi}

#' Convert degrees into radians
#' @noRd
deg2rad <- function(deg) {(deg * pi) / 180}

#' Convert dbh into basal area
#' @noRd
dbh2ba <- function(x) {pi/4*x*x}

#' Convert basal area into dbh
#' @noRd
ba2dbh <- function(x) {sqrt(x*4/pi)}

#' Keep value between lower and upper
#' @noRd
bracket <- function(x, lower, upper) {
  x[x >= lower & x <= upper]
}


#' Check if 3d points are equals at epsilon
#' @noRd
is_equal_3d <- function(x1, y1, z1,
                        x2, y2, z2,
                        epsilon) {

  abs(x1 - x2) < epsilon &
    abs(y1 - y2) < epsilon &
    abs(z1 - z2) < epsilon

}


#' Compute bottom azimut
#' 
#' @param north2x double - Angle from North to x axis clockwise. (in degrees)
#'    Default correspond to a Y axis oriented toward the North.
#' @param aspect double - Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
#'    northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
#'    
#' @export
#' @keywords internals
#' 
get_bottom_azimut <- function(aspect, 
                              north2x) {
  
  -aspect + north2x
  
}



#' Compute z coordinate of a point (x,y).
#'
#' @param x X-coordinate of the point
#' @param y Y-coordinate of the point
#' @param slope_rad Slope of the stand (in radians)
#' @param bottom_azimut_rad Azimuth of the vector orthogonal to the ground in the x,y system (in radians)
#' 
#' @importFrom data.table fifelse
#' 
#' @export
#' @keywords internals
#' 
get_z <- function(x, y,
                  slope_rad,
                  bottom_azimut_rad) {

  # Distance to origin
  d <- sqrt(x * x + y * y)

  # Azimut of the point (x,y,z) in the reference system
  azimut_rad <- data.table::fifelse(d == 0, 0,
                                    data.table::fifelse(y >= 0,
                                                        acos(x / d),
                                                        2 * pi - acos(x / d)))


  -d * cos(azimut_rad - bottom_azimut_rad) * tan(slope_rad)
}


#' Get distance between point 1 (x1, y1) and point 2 (x2, y2)
#' @noRd
get_distance <- function(x1, y1,
                         x2, y2) {

  sqrt( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))

}


#' Convert relative to absolute cell position
#'
#' @param x_target X-coordinate of the target cell (in meters)
#' @param y_target Y-coordinate of the target cell (in meters)
#' @param xid_rel Relative X-id from the target cell
#' @param yid_rel Relative Y-id from the target cell
#' @param cell_size Length of the side of a squared cell composing the stand (in meters)
#' @param n_cells integer - Number of cells of the side of a squared stand
#'
#' @return A list with:
#'  \itemize{
#'    \item{"xid_neighbour_eq": }{X-id of the neighbour cell (or equivalent within the grid if cell is outside the grid)}
#'    \item{"yid_neighbour_eq": }{Y-id of the neighbour cell (or equivalent within the grid if cell is outside the grid)}
#'    \item{"x_shift": }{Shift in X-coordinate from equivalent cell if the cell is outsie the grid (torus sytem)}
#'    \item{"y_shift": }{Shift in Y-coordinate from equivalent cell if the cell is outsie the grid (torus sytem)}
#'  }
#' @noRd
cells_rel2abs <- function(x_target, y_target,
                          xid_rel, yid_rel,
                          cell_size, n_cells) {

  # Transform target cell coordinates into id (be careful, Y-coordinate system is opposite direction to Y-grid system)
  xid_target <- x_target %/% cell_size
  yid_target <- n_cells - y_target %/% cell_size - 1


  # Get absolute id of neighbour cells
  xid_neighbour <- xid_target + xid_rel
  yid_neighbour <- yid_target + yid_rel


  # Search if cells is outside the plot
  outside <-
    xid_neighbour < 0 | xid_neighbour >= n_cells |
    yid_neighbour < 0 | yid_neighbour >= n_cells


  # Get equivalent in the grid cell using modulo (within a torus system)
  xid_neighbour_eq <- xid_neighbour %% n_cells
  yid_neighbour_eq <- yid_neighbour %% n_cells


  # Compute shift in coordinates to apply on trees within the equivalent cell (torus system)
  # Shift is a multiple of plot_size: how many plot size we have to shift tree coordinates times plot size, to apply torus system
  # If we have a neighbour cell in negative coordinates (id_neighbour_eq is always between [0, n_cells[ but negative id_neighbour)
  # ==> Thus, we want to substract coordinates of equivalent cells in order to have neighbour cell within negative coordinates
  # Be careful, y-coordinates system is opposite direction of y-grid system
  plot_size <- cell_size * n_cells

  x_shift <- (xid_neighbour %/% n_cells) * plot_size # Negative x_id ==> negative shift
  y_shift <- - (yid_neighbour %/% n_cells) * plot_size # Negative y id ==> positive shift


  # Output list
  list(
    xid_neighbour_eq, yid_neighbour_eq,
    x_shift, y_shift, outside
  )
}
