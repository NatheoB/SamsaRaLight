#' Print a \code{sl_stand} object
#'
#' Prints a compact, human-readable, one-line description of a
#' \code{sl_stand} object, summarizing the stand size, number of trees,
#' number of sensors, and grid geometry.
#'
#' This method is automatically called when a \code{sl_stand} object
#' is printed in the console.
#'
#' @param x An object of class \code{sl_stand}.
#' @param ... Unused. Included for S3 method compatibility.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print sl_stand
#' @export
#'
print.sl_stand <- function(x, ...) {
  
  stopifnot(inherits(x, "sl_stand"))
  
  trees <- x$trees
  sensors <- x$sensors
  geom  <- x$geometry
  tr    <- x$transform
  
  # ---- Tree count ----
  total_trees <- nrow(trees)
  
  # ---- Stand area ----
  total_area <- tr$new_area_ha
  
  # ---- Sensors count ----
  total_sensors <- ifelse(is.null(nrow(sensors)), 0, nrow(sensors))
  
  # ---- Geometry ----
  nx <- geom$n_cells_x
  ny <- geom$n_cells_y
  cs <- geom$cell_size
  
  # ---- Build label ----
  label <- "SamsaRaLight"
  
  
  label <- paste0(
    label,
    " stand of ",
    format(total_area, digits = 3),
    " ha with ",
    total_trees,
    " trees and ",
    total_sensors,
    " sensors"
  )
  
  label <- paste0(
    label,
    " (",
    nx, " x ", ny, " cells, ",
    cs, " m)"
  )
  
  cat(label, "\n")
  
  invisible(x)
}
