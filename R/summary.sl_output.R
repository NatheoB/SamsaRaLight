#' Summary of a SamsaRaLight simulation
#'
#' @param object an object of class \code{sl_output}
#' @param ... unused
#'
#' @export
#' @method summary sl_output
#'
summary.sl_output <- function(object, ...) {
  
  stopifnot(inherits(object, "sl_output"))
  
  cat("\n")
  cat("SamsaRaLight simulation summary\n")
  cat("================================\n\n")
  
  light  <- object$output$light
  params <- object$input$params
  
  # ==============================
  # TREE SUMMARY
  # ==============================
  cat("Trees (crown interception)\n")
  cat("---------------------------\n")
  
  tree_vars <- c("epot", "e", "lci")
  tree_vars <- tree_vars[tree_vars %in% names(light$trees)]
  
  if (length(tree_vars)) {
    print(summary(light$trees[, tree_vars, drop = FALSE]))
  } else {
    cat("No tree energy variables available\n")
  }
  cat("\n")
  
  
  # ==============================
  # CELL SUMMARY
  # ==============================
  cat("Cells (ground light)\n")
  cat("-------------------\n")
  
  cell_vars <- c("e", "pacl", "punobs")
  cell_vars <- cell_vars[cell_vars %in% names(light$cells)]
  
  if (length(cell_vars)) {
    print(summary(light$cells[, cell_vars, drop = FALSE]))
  } else {
    cat("No cell energy variables available\n")
  }
  cat("\n")
  
  
  # ==============================
  # SENSOR SUMMARY
  # ==============================
  cat("Sensors\n")
  cat("-------\n")
  
  sensor_vars <- c("e", "pacl", "punobs")
  sensor_vars <- sensor_vars[sensor_vars %in% names(light$sensors)]
  
  if (length(sensor_vars) && nrow(light$sensors) > 0) {
    print(summary(light$sensors[, sensor_vars, drop = FALSE]))
  } else {
    cat("No sensor energy variables available\n")
  }
  cat("\n")
  
  
  invisible(object)
}
