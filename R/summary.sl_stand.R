#' Summary of a SamsaRaLight stand
#'
#' @param object A \code{sl_stand} object
#' @param ... Unused
#'
#' @return Invisibly returns a list of summary tables
#' 
#' @export
#' @method summary sl_stand
#' 
summary.sl_stand <- function(object, ...) {
  
  stopifnot(inherits(object, "sl_stand"))
  
  trees <- object$trees
  sensors <- object$sensors
  geom  <- object$geometry
  trans <- object$transform
  
  ## -------------------------------
  ## Split trees
  ## -------------------------------
  trees_core  <- trees[!trees$added_to_fill, ]
  trees_all   <- trees
  
  ## -------------------------------
  ## Stand area
  ## -------------------------------
  area_core_ha <- trans$core_area_ha
  area_all_ha  <- trans$new_area_ha
  
  ## -------------------------------
  ## Forestry metrics
  ## -------------------------------
  stand_metrics <- function(tr, area_ha) {
    if (nrow(tr) == 0) {
      return(c(
        n = 0,
        density = NA,
        dg_cm = NA,
        ba_m2ha = NA
      ))
    }
    
    dbh <- tr$dbh_cm
    
    ba <- sum(pi * (dbh / 200)^2)      # m2
    dg <- sqrt(mean(dbh^2))            # quadratic mean diameter
    
    c(
      n = nrow(tr),
      density = nrow(tr) / area_ha,
      dg_cm = dg,
      ba_m2ha = ba / area_ha
    )
  }
  
  core_stats <- stand_metrics(trees_core, area_core_ha)
  all_stats  <- stand_metrics(trees_all,  area_all_ha)
  
  ## -------------------------------
  ## Geometry
  ## -------------------------------
  geom_stats <- list(
    cell_size = geom$cell_size,
    n_cells_x = geom$n_cells_x,
    n_cells_y = geom$n_cells_y,
    n_cells   = geom$n_cells_x * geom$n_cells_y,
    slope     = geom$slope,
    aspect    = geom$aspect,
    north2x   = geom$north2x
  )
  
  ## -------------------------------
  ## Sensors
  ## -------------------------------
  n_sensors <- ifelse(is.null(nrow(sensors)), 0, nrow(sensors))
  
  ## -------------------------------
  ## Light interception model
  ## -------------------------------
  lad_ok   <- all(!is.na(trees$crown_lad))
  open_ok  <- all(!is.na(trees$crown_openness))
  
  model <- ""
  if (lad_ok) {
    model <- "  - Turbid medium (crown_lad)"
  } 
  if (open_ok) {
    str_porous <- "  - Porous envelope (crown_openness)"
    model <- ifelse(lad_ok, paste(model, str_porous, sep = "\n"), str_porous) 
  }
  
  if (model == "") {
    model <- "  !!! No valid interception model (missing crown properties)"
  }
  
  ## -------------------------------
  ## Print
  ## -------------------------------
  cat("\n")
  cat("SamsaRaLight stand summary\n")
  cat("================================\n\n")
  
  cat("\nInventory (core polygon):\n")
  cat(sprintf("  Area              : %.2f ha\n", area_core_ha))
  cat(sprintf("  Trees             : %d\n", core_stats["n"]))
  cat(sprintf("  Density           : %.1f trees/ha\n", core_stats["density"]))
  cat(sprintf("  Basal area        : %.2f m2/ha\n", core_stats["ba_m2ha"]))
  cat(sprintf("  Quadratic mean DBH: %.1f cm\n", core_stats["dg_cm"]))
  
  cat("\nSimulation stand (core + filled buffer):\n")
  cat(sprintf("  Area              : %.2f ha\n", area_all_ha))
  cat(sprintf("  Trees             : %d\n", all_stats["n"]))
  cat(sprintf("  Density           : %.1f trees/ha\n", all_stats["density"]))
  cat(sprintf("  Basal area        : %.2f m2/ha\n", all_stats["ba_m2ha"]))
  cat(sprintf("  Quadratic mean DBH: %.1f cm\n", all_stats["dg_cm"]))
  
  cat("\nStand geometry:\n")
  cat(sprintf("  Grid              : %d x %d (%d cells)\n",
              geom$n_cells_x, geom$n_cells_y, geom_stats$n_cells))
  cat(sprintf("  Cell size         : %.2f m\n", geom$cell_size))
  cat(sprintf("  Slope             : %.2f deg\n", geom$slope))
  cat(sprintf("  Aspect            : %.2f deg\n", geom$aspect))
  cat(sprintf("  North to X-axis   : %.2f deg\n", geom$north2x))
  
  cat("\nNumber of sensors: ")
  cat(sprintf("%s\n", n_sensors))
  
  cat("\nAvailable light interception models:\n")
  cat(sprintf("%s\n", model))
  
  cat("\n")
  
  invisible(list(
    core = core_stats,
    full = all_stats,
    geometry = geom_stats,
    model = model
  ))
}
