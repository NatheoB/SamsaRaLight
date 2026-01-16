#' Print a \code{sl_output} object
#'
#' Prints a concise one-line summary of a SamsaRaLight simulation result,
#' reporting the number of cells, trees, and sensors, and whether detailed
#' outputs are available.
#'
#' @param x An object of class \code{sl_output}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return The object \code{x}, invisibly.
#'
#' @export
#' @method print sl_output
print.sl_output <- function(x, ...) {
  
  stopifnot(inherits(x, "sl_output"))
  
  stand <- x$input$sl_stand
  
  n_cells   <- if (!is.null(stand$cells))   nrow(stand$cells)   else 0
  n_trees   <- if (!is.null(stand$trees))   nrow(stand$trees)   else 0
  n_sensors <- if (!is.null(stand$sensors)) nrow(stand$sensors) else 0
  
  # Detailed output = rays + interceptions
  has_detailed <- !is.null(x$output$monthly_rays) &&
    !is.null(x$output$interceptions)
  
  cat(
    "SamsaRaLight output with",
    n_cells, "cells,",
    n_trees, "trees and",
    n_sensors, "sensors",
    if (has_detailed) "(detailed output)\n" else "(no detailed output)\n"
  )
  
  invisible(x)
}
