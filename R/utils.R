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
#' @return a bottom azimuth numeric value
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
#' @return a z coordinate numeric value
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

#' Rotate 2D vectors by an angle (counter-clockwise)
#'
#' Internal utility to rotate 2D coordinates by a given angle (in radians),
#' using a standard rotation matrix.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates.
#' @param theta Rotation angle in radians (counter-clockwise).
#'
#' @return A data.frame with rotated coordinates:
#' \describe{
#'   \item{x}{Rotated x coordinates}
#'   \item{y}{Rotated y coordinates}
#' }
#'
#' @details
#' The rotation follows:
#' \deqn{
#' x' = x \cos(\theta) - y \sin(\theta)
#' }
#' \deqn{
#' y' = x \sin(\theta) + y \cos(\theta)
#' }
#'
#' @keywords internal
rotate_vec_ccw <- function(x, y, theta) {
  data.frame(
    x = x * cos(theta) - y * sin(theta),
    y = x * sin(theta) + y * cos(theta)
  )
}


#' Check the arguments common to the SamsaRaLight runs
#'
#' Validates the stand, the monthly radiations, the logical and numeric scalar
#' arguments, the simulated period and the number of OpenMP threads
#' (used by \link{run_sl_advanced} and \link{precompute_sl_geometry}).
#'
#' @param sl_stand An object of class \code{"sl_stand"}.
#' @param monthly_radiations A data frame of monthly radiations.
#' @param logical_args Named list of the arguments that must be a single TRUE or FALSE.
#' @param numeric_args Named list of the arguments that must be a single numeric value.
#' @param start_day,end_day First and last day of the simulated period.
#' @param n_threads Number of OpenMP threads, or \code{NULL}.
#'
#' @return Invisibly \code{TRUE}, or an error.
#'
#' @noRd
check_sl_run_args <- function(sl_stand,
                              monthly_radiations,
                              logical_args,
                              numeric_args,
                              start_day,
                              end_day,
                              n_threads) {

  # sl_stand
  validate_sl_stand(sl_stand)

  # monthly_radiations
  check_monthly_radiations(monthly_radiations, verbose = FALSE)

  # Logical arguments
  for (arg_name in names(logical_args)) {
    arg_val <- logical_args[[arg_name]]
    if (!is.logical(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single TRUE or FALSE.", arg_name), call. = FALSE)
    }
  }

  # Numeric scalar arguments
  for (arg_name in names(numeric_args)) {
    arg_val <- numeric_args[[arg_name]]
    if (!is.numeric(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single numeric value.", arg_name), call. = FALSE)
    }
  }

  # Days
  if (!is.numeric(start_day) || length(start_day) != 1 || start_day < 1 || start_day > 365) {
    stop("`start_day` must be a single number between 1 and 365.", call. = FALSE)
  }
  if (!is.numeric(end_day) || length(end_day) != 1 || end_day < 1 || end_day > 365) {
    stop("`end_day` must be a single number between 1 and 365.", call. = FALSE)
  }
  if (end_day < start_day) {
    stop("`end_day` cannot be smaller than `start_day`.", call. = FALSE)
  }

  # Parallel threads
  if (!is.null(n_threads)) {
    if (!is.numeric(n_threads) || length(n_threads) != 1 || is.na(n_threads)) {
      stop("`n_threads` must be a single integer or NULL.", call. = FALSE)
    }
    if (n_threads < 1 || n_threads %% 1 != 0) {
      stop("`n_threads` must be a positive integer.", call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' Set the parallel computation of the SamsaRaLight C++ code
#'
#' Forces BLAS and R-level OpenMP to a single thread (to avoid competing
#' parallel execution), then sets the OpenMP threads used by the C++ code.
#'
#' @param parallel_mode Logical. Parallelise with OpenMP?
#' @param n_threads Number of OpenMP threads, or \code{NULL} (automatic).
#' @param verbose Logical. Print the OpenMP status?
#'
#' @return Invisibly \code{NULL}.
#'
#' @noRd
set_sl_parallel <- function(parallel_mode, n_threads, verbose) {

  ## Prevent BLAS and OpenMP fighting ----
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)

  ## Set parallel ----
  if (is.null(n_threads)) n_threads <- -1L
  sl_set_openmp(parallel_mode = parallel_mode,
                num_threads = as.integer(n_threads),
                verbose = verbose)

  if (verbose) sl_print_openmp_status()

  invisible(NULL)
}


#' Warn about the interceptions ignored by the C++ code
#'
#' Interceptions with an odd number of interception points that are not a
#' known particular case are ignored, as in Capsis, and reported to the user
#' in a single warning.
#'
#' @param n_unexpected Named numeric vector of the number of ignored
#'   interceptions (\code{ellipsoid}, \code{paraboloid}, \code{trunk}).
#'
#' @return Invisibly \code{NULL}.
#'
#' @noRd
warn_sl_unexpected_interceptions <- function(n_unexpected) {

  if (sum(n_unexpected) > 0) {
    warning(sprintf(
      paste0("%s ray/tree interception(s) with an unexpected number of interception points ",
             "were ignored, as in Capsis (ellipsoids: %s, paraboloids: %s, trunks: %s). ",
             "This may come from rays crossing crown or trunk edges."),
      sum(n_unexpected), n_unexpected[["ellipsoid"]], n_unexpected[["paraboloid"]],
      n_unexpected[["trunk"]]),
      call. = FALSE)
  }

  invisible(NULL)
}
