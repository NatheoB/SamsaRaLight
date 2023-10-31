#' Compute ray interception with an ellipsoidal crown
#'
#' @description Evaluates if the given ray intercepts this crown, previously relocated
#'  with the given shifts.
#'
#' @param cos_elevation double - Cosinus ofelevation angle of the ray (in radians)
#' @param sin_elevation double - Sinus ofelevation angle of the ray (in radians)
#' @param cos_azimut double - Cosinus of azimut angle of the ray (in radians, trigonometric 0 and counter clockwise)
#' @param sin_azimut double - Sinus of azimut angle of the ray (in radians, trigonometric 0 and counter clockwise)
#' @param x double - X-component of paraboloid center (in meters)
#' @param y double - Y-component of paraboloid center (in meters)
#' @param z double - Paraboloid base height (in meters)
#' @param x_shift double - Shift of paraboloid center on X-axis (in meters)
#' @param y_shift double - Shift of paraboloid center on Y-axis (in meters)
#' @param z_shift double - Shift of paraboloid center on Z-axis (in meters)
#' @param a double - Parameter a of paraboloid = semi-axis of basal ellipse along X axis
#' @param b double - Parameter b of paraboloid = semi-axis of basal ellipse along Y axis
#' @param h double - Parameter h of paraboloid = crown height, along Z axis
#'
#' @return Returns the interception path length and a distance to
#'  the origin (i.e. the target cell center) (NA is no interception)
#'
#' @importFrom matrixStats rowMaxs rowMins
#'
#' @export
sl_intercept_crown_paraboloid <- function(cos_elevation, sin_elevation,
                                          cos_azimut, sin_azimut,
                                          x, y, z,
                                          x_shift, y_shift, z_shift,
                                          a, b, h) {


  ### GET SHIFTED POSITION OF THE PARABOLOID ----
  x <- x + x_shift
  y <- y + y_shift
  z <- z + z_shift


  ### FIND SOLUTION OF THE QUADRATIC EQUATION (A*x*x + B*x + C = 0)
  # Equation giving distance between ray intersection with tree crown and target cell center

  # Compute a, b and c coefficients
  coef_A <- cos_elevation * cos_elevation * cos_azimut * cos_azimut / (a * a) +
    cos_elevation * cos_elevation * sin_azimut * sin_azimut / (b * b)

  coef_B <- -2 * x * cos_elevation * cos_azimut / (a * a) -
    2 * y * sin_azimut * cos_elevation / (b * b) +
    sin_elevation / h

  coef_C <- x * x / (a * a) + y * y / (b * b) - (z + h) / h


  # Find matrix of positive solutions (distance to target cell), if not 2 solutions, set to NA
  # Because negative  = no interception and null = 1 interception = do not consider tangent rays with crown
  dist_mat <- solve_quadratic_equation_onlypos(coef_A, coef_B, coef_C)

  dist_mat <- cbind(dist_mat,
                    z / sin_elevation,
                    0)

  ### COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS ----

  # Interception points of the two roots
  v1 <- get_intercept_coords(dist_mat[,1],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  v2 <- get_intercept_coords(dist_mat[,2],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)

  # Interception with the base plane
  vbase <- get_intercept_coords(dist_mat[,3],
                                cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  vbase[,3] <- z # Correct for approximations

  ### FIND PATH LENGTH AND DISTANCE TO CELL CENTER ----

  # Get the limits of the box between shifted paraboloid center and crown limit
  x_bbox_max <- x + a
  y_bbox_max <- y + b
  z_bbox_max <- z + h

  x_bbox_min <- x - a
  y_bbox_min <- y - b
  z_bbox_min <- z

  # Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)
  epsilon <- 1e-10 # For rounding errors

  is_v1 <- v1[,3]>=0 &
    in_bbox(v1[,1], v1[,2], v1[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  is_v2 <- v2[,3]>=0 &
    in_bbox(v2[,1], v2[,2], v2[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  is_vbase <- vbase[,3]>=0 &
    in_bbox(vbase[,1], vbase[,2], vbase[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_paraboloid(vbase[,1], vbase[,2], vbase[,3],
                  x, y, z, a, b, h)

  is_target <-
    in_bbox(0, 0, 0,
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_paraboloid(0, 0, 0,
                  x, y, z, a, b, h)


  # Find number of solutions
  is_sols <- cbind(is_v1, is_v2, is_vbase, is_target)

  n_sols <- rowSums(is_sols, na.rm = T)
  if (!all(n_sols %in% c(0, 2))) stop("Uncorrect number of interceptions")

  # Return pathlength and distance of mid path from target if ray has crossed the ellipsoid
  is_sols[!is_sols] <- NA
  dist_mat <- is_sols * dist_mat

  length <- matrixStats::rowMaxs(dist_mat, na.rm = TRUE) - matrixStats::rowMins(dist_mat, na.rm = TRUE)
  # length[length == -Inf] <- NA

  distance <- rowMeans(dist_mat, na.rm = T)
  # distance[is.nan(distance)] <- NA

  out <- list(length, distance)

  return(out)
}


#' Compute ray interception with an ellipsoidal crown
#'
#' @description Evaluates if the given ray intercepts this crown, previously relocated
#'  with the given shifts.
#'
#' @param cos_elevation double - Cosinus ofelevation angle of the ray (in radians)
#' @param sin_elevation double - Sinus ofelevation angle of the ray (in radians)
#' @param cos_azimut double - Cosinus of azimut angle of the ray (in radians, trigonometric 0 and counter clockwise)
#' @param sin_azimut double - Sinus of azimut angle of the ray (in radians, trigonometric 0 and counter clockwise)
#' @param x double - X-component of ellipsoid center (in meters)
#' @param y double - Y-component of ellipsoid center (in meters)
#' @param z double - Z-component of ellipsoid center (in meters)
#' @param x_shift double - Shift of ellipsoid center on X-axis (in meters)
#' @param y_shift double - Shift of ellipsoid center on Y-axis (in meters)
#' @param z_shift double - Shift of ellipsoid center on Z-axis (in meters)
#' @param a double - Parameter a of ellipsoid semi-principal axes (X-axis, semi-major)
#' @param b double - Parameter b of ellipsoid semi-principal axes (Y-axis, semi-minor)
#' @param c double - Parameter c of ellipsoid semi-principal axes (Z-axis, height)
#'
#' @return Returns the interception path length and a distance to
#'  the origin (i.e. the target cell center) (NA is no interception)
#'
#' @export
sl_intercept_crown_ellipsoid <- function(cos_elevation, sin_elevation,
                                         cos_azimut, sin_azimut,
                                         x, y, z,
                                         x_shift, y_shift, z_shift,
                                         a, b, c) {

  # Check for length of rays properties and ellipsoids

  ### GET SHIFTED POSITION OF THE ELLIPSOID ----
  x <- x + x_shift
  y <- y + y_shift
  z <- z + z_shift


  ### FIND SOLUTION OF THE QUADRATIC EQUATION (A*x*x + B*x + C = 0)
  # Equation giving distance between ray intersection with tree crown and target cell center

  # Compute a, b and c coefficients
  coef_A <- cos_elevation * cos_elevation * cos_azimut * cos_azimut / (a * a) +
    cos_elevation * cos_elevation * sin_azimut * sin_azimut / (b * b) +
    sin_elevation * sin_elevation / (c * c)

  coef_B <- -2 * x * cos_elevation * cos_azimut / (a * a) -
    2 * y * sin_azimut * cos_elevation / (b * b) -
    2 * z * sin_elevation / (c * c)

  coef_C <- x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - 1


  # Find matrix of positive solutions (distance to target cell)
  # Because negative  = no interception and null = 1 interception = do not consider tangent rays with crown
  dist_mat <- solve_quadratic_equation_onlypos(coef_A, coef_B, coef_C)

  ### COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS ----

  # Interception points of the two roots
  v1 <- get_intercept_coords(dist_mat[,1],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  v2 <- get_intercept_coords(dist_mat[,2],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)


  ### FIND PATH LENGTH AND DISTANCE TO CELL CENTER ----

  # Get the limits of the box between shifted ellispoid center and crown limit
  x_bbox_max <- x + a
  y_bbox_max <- y + b
  z_bbox_max <- z + c

  x_bbox_min <- x - a
  y_bbox_min <- y - b
  z_bbox_min <- z - c

  # Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)
  epsilon <- 1e-10 # For rounding errors

  is_v1 <- v1[,3]>=0 &
    in_bbox(v1[,1], v1[,2], v1[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  is_v2 <- v2[,3]>=0 &
    in_bbox(v2[,1], v2[,2], v2[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  # Find number of solutions
  is_sols <- cbind(is_v1, is_v2)

  n_sols <- rowSums(is_sols, na.rm = T)
  if (!all(n_sols %in% c(0, 2))) stop("Uncorrect number of interceptions")

  # Return pathlength and distance of mid path from target if ray has crossed the ellipsoid
  dist_mat <- is_sols * dist_mat

  length <- abs(dist_mat[,1] - dist_mat[,2])
  distance <- (dist_mat[,1] + dist_mat[,2]) / 2

  res <- list(length, distance)

  return(res)
}


#' Compute ray interception with a 8th of ellipsoidal crown
#'
#' @description Evaluates if the given ray intercepts this crown, previously relocated
#'  with the given shifts.
#'
#' @param elevation double - Elevation angle of the ray (in radians)
#' @param azimut double - azimut angle of the ray (in radians, trigonometric 0 and counter clockwise)
#' @param x double - X-component of ellipsoid center (in meters)
#' @param y double - Y-component of ellipsoid center (in meters)
#' @param z double - Z-component of ellipsoid center (in meters)
#' @param x_shift double - Shift of ellipsoid center on X-axis (in meters)
#' @param y_shift double - Shift of ellipsoid center on Y-axis (in meters)
#' @param z_shift double - Shift of ellipsoid center on Z-axis (in meters)
#' @param a double - Parameter a of ellipsoid semi-principal axes (X-axis, semi-major)
#' @param b double - Parameter b of ellipsoid semi-principal axes (Y-axis, semi-minor)
#' @param c double - Parameter c of ellipsoid semi-principal axes (Z-axis, height)
#'
#' @return Returns the interception path length and a distance to
#'  the origin (i.e. the target cell center) (NA is no interception)
#'
#' @export
sl_intercept_crown_8thellipsoidal <- function(elevation, azimut,
                                              x, y, z,
                                              x_shift, y_shift, z_shift,
                                              a, b, c) {

  # Check for length of rays properties and ellipsoids


  ### GET SHIFTED POSITION OF THE ELLIPSOID ----
  x <- x + x_shift
  y <- y + y_shift
  z <- z + z_shift


  ### FIND SOLUTION OF THE QUADRATIC EQUATION (A*x*x + B*x + C = 0)
  # Equation giving distance between ray intersection with tree crown and target cell center

    # Compute a, b and c coefficients
  cos_elevation <- cos(elevation)
  sin_elevation <- sin(elevation)
  cos_azimut <- cos(azimut)
  sin_azimut <- sin(azimut)

  coef_A <- cos_elevation * cos_elevation * cos_azimut * cos_azimut / (a * a) +
    cos_elevation * cos_elevation * sin_azimut * sin_azimut / (b * b) +
    sin_elevation * sin_elevation / (c * c)

  coef_B <- -2 * x * cos_elevation * cos_azimut / (a * a) -
    2 * y * sin_azimut * cos_elevation / (b * b) -
    2 * z * sin_elevation / (c * c)

  coef_C <- x * x / (a * a) + y * y / (b * b) + z * z / (c * c) - 1


    # Find matrix of positive solutions (lengths)
    # Because negative  = no interception and null = 1 interception = do not consider tangent rays with crown
  dist_mat <- solve_quadratic_equation_onlypos(coef_A, coef_B, coef_C)

    # Add lengths of intercepted points with plan x, y, z and target cell (i.e. always 0 because origin)
  dist_mat <- cbind(dist_mat,
                    x / (cos_elevation * cos_azimut),
                    y / (cos_elevation * sin_azimut),
                    z / sin_elevation,
                    0)


  ### COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS ----

    # Interception points of the two roots
  v1 <- get_intercept_coords(dist_mat[,1],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  v2 <- get_intercept_coords(dist_mat[,2],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)

  # Interception with plane x
  vx <- get_intercept_coords(dist_mat[,3],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)

  vx[,1] <- x #for rounding errors

    # Interception with plane y
  vy <- get_intercept_coords(dist_mat[,4],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  vy[,2] <- y #for rounding errors

  # Interception with plane z
  vz <- get_intercept_coords(dist_mat[,5],
                             cos_elevation, sin_elevation, cos_azimut, sin_azimut)
  vz[,3] <- z #for rounding errors


  ### FIND PATH LENGTH AND DISTANCE TO CELL CENTER ----

    # Top limits of ellipsoid
  xa <- x + a
  yb <- y + b
  zc <- z + c

    # Get the limits of the box between shifted ellispoid center and crown limit
  x_bbox_min <- pmin(x, xa)
  y_bbox_min <- pmin(y, yb)
  z_bbox_min <- pmin(z, zc)
  x_bbox_max <- pmax(x, xa)
  y_bbox_max <- pmax(y, yb)
  z_bbox_max <- pmax(z, zc)


    # Find interception points within the crown
  epsilon <- 1e-10 # For rounding errors

  is_v1 <- v1[,3]>=0 &
    in_bbox(v1[,1], v1[,2], v1[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  is_v2 <- v2[,3]>=0 &
    in_bbox(v2[,1], v2[,2], v2[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max)

  is_vx <- vx[,3]>0 &
    in_bbox(vx[,1], vx[,2], vx[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_ellipsoid(vx[,1], vx[,2], vx[,3],
                 x, y, z, a, b, c)

  is_vy <- vy[,3]>0 &
    in_bbox(vy[,1], vy[,2], vy[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_ellipsoid(vy[,1], vy[,2], vy[,3],
                 x, y, z, a, b, c) &
    !is_equal_3d(vx[,1], vx[,2], vx[,3],
                 vy[,1], vy[,2], vy[,3],
                 epsilon)

  is_vz <- vz[,3]>0 &
    in_bbox(vz[,1], vz[,2], vz[,3],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_ellipsoid(vz[,1], vz[,2], vz[,3],
                 x, y, z, a, b, c) &
    !is_equal_3d(vz[,1], vz[,2], vz[,3],
                 vx[,1], vx[,2], vx[,3],
                 epsilon) &
    !is_equal_3d(vz[,1], vz[,2], vz[,3],
                 vy[,1], vy[,2], vy[,3],
                 epsilon)

  is_target <-
    in_bbox(dist_mat[,6], dist_mat[,6], dist_mat[,6],
            x_bbox_min, y_bbox_min, z_bbox_min,
            x_bbox_max, y_bbox_max, z_bbox_max) &
    in_ellipsoid(dist_mat[,6], dist_mat[,6], dist_mat[,6],
                 x, y, z, a, b, c)


    # Find number of solutions
  is_sols <- cbind(
    is_v1, is_v2,
    is_vx, is_vy, is_vz,
    is_target
  )
  is_sols[!is_sols] <- NA

  n_sols <- rowSums(is_sols, na.rm = T)
  if (!all(n_sols %in% c(0, 2))) stop("Uncorrect number of interceptions")

  # Return pathlength and distance of mid path from target if ray has crossed the ellipsoid
  dist_mat <- is_sols * dist_mat

  length <- matrixStats::rowMaxs(dist_mat, na.rm = TRUE) - matrixStats::rowMins(dist_mat, na.rm = TRUE)
  # length[length == -Inf] <- NA

  distance <- rowMeans(dist_mat, na.rm = T)
  # distance[is.nan(distance)] <- NA

  out <- list(length, distance)

  return(out)
}


#' Solve a quadratic equation (ax^2 + bx + c) and only compute positive delta
#' @return Return a matrix with the 2 solutions as columns (NULL values for negative delta)
#' @noRd
#' @export
solve_quadratic_equation_onlypos <- function(a, b, c) {

  # Compute discriminant
  delta <- b * b - 4 * a * c

  # Init output
  sol_mat <- matrix(NA, nrow = length(a), ncol = 2)

  # Compute sqrt of discriminant if positive
  delta_pos <- delta > 0

  delta_sqrt <- sqrt(delta[delta_pos])
  a <- a[delta_pos]
  b <- b[delta_pos]

  # Get two solutions if delta pos
  sol_mat[delta_pos,1] <- (- b + delta_sqrt) / (2 * a)
  sol_mat[delta_pos,2] <- (- b - delta_sqrt) / (2 * a)

  return(sol_mat)
}


#' Compute the coordinates of the interception point from the distance to
#' the target cell and beam direction.
#'
#' @param length - Distance from intersection point and target center cell (in meters)
#' @param cos_elevation - Cosinus of elevation angle of the beam (in radians)
#' @param sin_elevation - Sinus of elevation angle of the beam (in radians)
#' @param cos_azimut - Cosinus of azimut angle of the beam (in radians)
#' @param sin_azimut - Sinus of azimut angle of the beam (in radians)
#'
#' @return a vector containing the coordinates of the interception point and
#' its distance to the target cell
#'
#' @export
get_intercept_coords <- function(length,
                                 cos_elevation, sin_elevation,
                                 cos_azimut, sin_azimut) {

    x <- length * cos_elevation * cos_azimut
    y <- length * cos_elevation * sin_azimut
    z <- length * sin_elevation

    matrix(data = c(x, y, z, length), ncol = 4, byrow = FALSE)
  }


#' Check if a point is contained within a box at epsilon
#'
#' @param x x-coordinate of the point (double)
#' @param y y-coordinate of the point (double)
#' @param z z-coordinate of the point (double)
#' @param x_bbox_min Minimum x-coordinate of the box (double)
#' @param y_bbox_min Minimum y-coordinate of the box (double)
#' @param z_bbox_min Minimum z-coordinate of the box (double)
#' @param x_bbox_max Maximum x-coordinate of the box (double)
#' @param y_bbox_max Maximum y-coordinate of the box (double)
#' @param z_bbox_max Maximum z-coordinate of the box (double)
#' @noRd
#' @export
in_bbox <- function(x, y, z,
                    x_bbox_min, y_bbox_min, z_bbox_min,
                    x_bbox_max, y_bbox_max, z_bbox_max) {

  x >= x_bbox_min & x <= x_bbox_max &
    y >= y_bbox_min & y <= y_bbox_max &
    z >= z_bbox_min & z <= z_bbox_max

}


#' Check if a point is contained within a box at epsilon
#'
#' @param x x-coordinate of the point (double)
#' @param y y-coordinate of the point (double)
#' @param z z-coordinate of the point (double)
#' @param x_bbox_min Minimum x-coordinate of the box (double)
#' @param y_bbox_min Minimum y-coordinate of the box (double)
#' @param z_bbox_min Minimum z-coordinate of the box (double)
#' @param x_bbox_max Maximum x-coordinate of the box (double)
#' @param y_bbox_max Maximum y-coordinate of the box (double)
#' @param z_bbox_max Maximum z-coordinate of the box (double)
#' @param epsilon Accepted error for rounding approximation (double)
#' @noRd
#' @export
in_bbox_epsilon <- function(x, y, z,
                            x_bbox_min, y_bbox_min, z_bbox_min,
                            x_bbox_max, y_bbox_max, z_bbox_max,
                            epsilon) {

  x_bbox_min - x <= epsilon & x - x_bbox_max <= epsilon &
    y_bbox_min - y <= epsilon & y - y_bbox_max <= epsilon &
    z_bbox_min - z <= epsilon & z - z_bbox_max <= epsilon

}

#' Check if a point is contained within a paraboloid 3D shape
#'
#' @param x x-coordinate of the point (double)
#' @param y y-coordinate of the point (double)
#' @param z z-coordinate of the point (double)
#' @param x_par double - X-component of paraboloid center (in meters)
#' @param y_par double - Y-component of paraboloid center (in meters)
#' @param z_par double - Z-component of paraboloid center (in meters)
#' @param a double - Parameter a of paraboloid = semi-axis of basal ellipse along X axis
#' @param b double - Parameter b of paraboloid = semi-axis of basal ellipse along Y axis
#' @param h double - Parameter h of paraboloid = crown height, along Z axis
#' @noRd
#' @export
in_paraboloid <- function(x, y, z,
                          x_par, y_par, z_par,
                          a, b, h) {

  val <- (x - x_par) * (x - x_par) / (a * a) +
    (y - y_par) * (y - y_par) / (b * b) -
    (z - z_par) / (h * h)

  return(val <= 1)
}


#' Check if a point is contained within an ellipsoid 3D shape
#'
#' @param x x-coordinate of the point (double)
#' @param y y-coordinate of the point (double)
#' @param z z-coordinate of the point (double)
#' @param x_el x-coordinate of the center of the ellipsoid (double)
#' @param y_el y-coordinate of the center of the ellipsoid (double)
#' @param z_el z-coordinate of the center of the ellipsoid (double)
#' @param a double - Parameter a of ellipsoid semi-principal axes (X-axis, semi-major)
#' @param b double - Parameter b of ellipsoid semi-principal axes (Y-axis, semi-minor)
#' @param c double - Parameter c of ellipsoid semi-principal axes (Z-axis, height)
#' @noRd
#' @export
in_ellipsoid <- function(x, y, z,
                         x_el, y_el, z_el,
                         a, b, c) {

  val <- (x - x_el) * (x - x_el) / (a * a) +
    (y - y_el) * (y - y_el) / (b * b) +
    (z - z_el) * (z - z_el) / (c * c)

  return(val < 1)
}
