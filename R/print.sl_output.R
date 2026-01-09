#' Print a \code{sl_output} object
#'
#' Prints a concise overview of a SamsaRaLight simulation result, including
#' stand size, grid geometry, model configuration, and which outputs
#' (trees, cells, sensors, rays, interceptions) are available.
#'
#' This method is designed to give a quick human-readable summary of what
#' was simulated and what was produced, without printing large data tables.
#'
#' @param x An object of class \code{sl_output}.
#' @param ... Further arguments passed to or from other methods (currently ignored).
#'
#' @return The object \code{x}, invisibly.
#'
#' @details
#' The printed output includes:
#' \itemize{
#'   \item Stand size (area, number of trees, cells, sensors)
#'   \item Grid geometry (number of cells and resolution)
#'   \item Radiation model parameters (latitude, time period, turbid medium, etc.)
#'   \item Which light and interception outputs are available
#' }
#'
#' Use \code{\link{summary.sl_output}} to obtain numerical summaries
#' of light, interception and energy variables, and \code{\link{plot.sl_output}}
#' for graphical inspection.
#'
#' @seealso
#' \code{\link{summary.sl_output}}, \code{\link{print.sl_stand}}
#'
#' @export
#' @method print sl_output
print.sl_output <- function(x, ...) {
  
  stopifnot(inherits(x, "sl_output"))
  
  stand  <- x$input$sl_stand
  params <- x$input$params
  geom   <- stand$geometry
  tr     <- stand$transform
  light  <- x$output$light
  
  # ---- Sizes ----
  n_trees   <- nrow(stand$trees)
  n_cells   <- nrow(stand$cells)
  n_sensors <- if (is.null(nrow(stand$sensors))) 0 else nrow(stand$sensors)
    
  
  area_ha <- tr$new_area_ha
  
  nx <- geom$n_cells_x
  ny <- geom$n_cells_y
  cs <- geom$cell_size
  
  # ---- Rays ----
  has_rays <- !is.null(x$output$monthly_rays)
  n_rays <- if (has_rays) nrow(x$output$monthly_rays$rays) else 0
  
  # ---- Output flags ----
  has_tree_light   <- !is.null(light$trees)
  has_cell_light   <- !is.null(light$cells)
  has_sensor_light <- !is.null(light$sensors)
  
  has_intercept <- !is.null(x$output$interceptions)
  
  # ---- Header ----
  cat("SamsaRaLight simulation\n\n")
  
  # ---- Stand ----
  cat("Stand\n")
  cat("  Area    :", format(area_ha, digits = 3), "ha\n")
  cat("  Trees   :", n_trees, "\n")
  cat("  Cells   :", n_cells, sprintf(" (%d x %d at %g m)", nx, ny, cs), "\n")
  cat("  Sensors :", n_sensors, "\n\n")
  
  # ---- Outputs ----
  cat("Outputs\n")
  cat("  Tree light   :", if (has_tree_light) paste0("yes (", nrow(light$trees), " trees)") else "no", "\n")
  cat("  Cell light   :", if (has_cell_light) paste0("yes (", nrow(light$cells), " cells)") else "no", "\n")
  cat("  Sensor light :", if (has_sensor_light) paste0("yes (", nrow(light$sensors), " sensors)") else "no", "\n")
  cat("  Monthly rays :", if (has_rays) paste0("yes (", n_rays, " rays)") else "no", "\n")
  cat("  Interception :", if (has_intercept) "yes (trees x sensors, trees x cells)" else "no", "\n")
  
  invisible(x)
}
