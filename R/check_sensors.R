#' Check the format and validity of a sensor position data.frame
#'
#' This function checks whether a sensor data.frame is correctly formatted
#' to be used as input for the ray-tracing model. It verifies the presence,
#' type, and validity of mandatory variables describing the position and
#' height of sensors within a forest stand.
#'
#' @param sensors A data.frame with one row per sensor and the following columns:
#' \itemize{
#'   \item{id_sensor}{Unique identifier of the sensor (numeric or character, no duplicates)}
#'   \item{x}{X position of the sensor (numeric)}
#'   \item{y}{Y position of the sensor (numeric)}
#'   \item{h_m}{Height above ground of the sensor (numeric, meters)}
#' }
#'
#' @param verbose Logical; if \code{TRUE}, informative messages are printed.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#'
#' @examples
#' \dontrun{
#' sensors <- data.frame(
#'   id_sensor = 1:3,
#'   x = c(10, 20, 30),
#'   y = c(5, 15, 25),
#'   h_m = c(1.5, 2.0, 1.8)
#' )
#'
#' check_sensors(sensors)
#' }
#'
#' @export
check_sensors <- function(sensors, verbose = TRUE) {
  
  ## ---- basic structure ------------------------------------------------------
  
  if (is.null(sensors)) {
    if (verbose) {
      message("Sensors successfully validated: NULL data.frame.")
    }
    return(invisible(TRUE))
  }
  
  if (!is.data.frame(sensors)) {
    stop("`sensors` must be a data.frame.", call. = FALSE)
  }
  
  if (nrow(sensors) == 0) {
    stop("`sensors` must contain at least one sensor or be NULL.", call. = FALSE)
  }
  
  ## ---- required columns -----------------------------------------------------
  required_cols <- c("id_sensor", "x", "y", "h_m")
  missing_cols <- setdiff(required_cols, names(sensors))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  ## ---- id uniqueness --------------------------------------------------------
  if (anyDuplicated(sensors$id_sensor)) {
    stop("`id` must contain unique values (no duplicates).", call. = FALSE)
  }
  
  ## ---- numeric checks -------------------------------------------------------
  numeric_cols <- c("x", "y", "h_m")
  non_numeric <- numeric_cols[!vapply(sensors[numeric_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "The following columns must be numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (any(sensors$h_m < 0, na.rm = TRUE)) {
    stop("Sensor height must be non-negative.", call. = FALSE)
  }
  
  ## ---- success --------------------------------------------------------------
  if (verbose) {
    message("Sensors successfully validated: structure and values are consistent.")
  }
  
  invisible(TRUE)
}
