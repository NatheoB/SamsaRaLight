#' Run SamsaraLight from Capsis samsaralightloader model
#'
#' @param capsis_folderpath Folderpath of the capsis source directory
#' @param inv_fp Filepath of the lilo inventory
#' @param meteo_fp Filepath of the lilo meteo file
#' @param export_dir Filepath of the export directory
#' @param run_cell Boolean, if TRUE, run analysis also for cell light interception
#' @param java_server_started Boolean, if False, no not start and clode java server 
#'  (to run multiple simulations one after the other without restarting java server)
#'
#' @noRd
sl_run_lilo <- function(capsis_folderpath,
                        inv_fp, meteo_fp, export_dir,
                        run_cell = TRUE,
                        java_server_started = FALSE) {
  
  # Connect to capsis
  if (!java_server_started) {
    setCapsisPath(path = capsis_folderpath)
    connectToCapsis()
  }
  
  # Creates and initializes the script.
  tree_dir <- createJavaObject("java.io.File", inv_fp)
  meteo_dir <- createJavaObject("java.io.File", meteo_fp)
  output_dir <- createJavaObject("java.io.File", export_dir)
  export_dir <- createJavaObject("java.io.File", export_dir)
  
  
  # Tree light interception
  export_tree <- createJavaObject("samsaralightloader.myscripts.LiloExportTreeLight", 
                                  tree_dir, meteo_dir, output_dir, export_dir)
  
  ## Make the simulation
  message("Doing analysis...")
  export_tree$doAnalysis()
  
  ## Export the results
  results_tree <- export_tree$getResults()
  indexes_tree <- as.integer(0:(results_tree$size()-1))
  items_tree <- results_tree$get(indexes_tree)
  capsis.output_tree <- tibble(id_tree = items_tree$id[indexes_tree+1],
                               e = items_tree$e[indexes_tree+1],
                               epot = items_tree$epot[indexes_tree+1]
  )

  
  capsis.output_cell <- NULL
  if (run_cell) { 
    
    # Cell interception
    export_cell <- createJavaObject("samsaralightloader.myscripts.LiloExportCellLight", 
                                    tree_dir, meteo_dir, output_dir, export_dir)
    
    ## Make the simulation
    export_cell$doAnalysis()
    
    ## Export the results
    results_cell<-export_cell$getResults()
    indexes_cell <- as.integer(0:(results_cell$size()-1))
    items_cell<-results_cell$get(indexes_cell)
    capsis.output_cell<-tibble(id = items_cell$id[indexes_cell+1],
                               x_center = items_cell$x[indexes_cell+1] + 2.5,
                               y_center = items_cell$y[indexes_cell+1] + 2.5,
                               e = items_cell$e[indexes_cell+1],
                               erel = items_cell$erel[indexes_cell+1] / 100)
  }
  
  # Unconnect from Capsis
  if (!java_server_started) {
    shutdownClient()
  }
  
  # Return list of outputs
  list("trees" = capsis.output_tree, "cells" = capsis.output_cell)
}



#' Run SamsaraLight from Capsis Samsara model
#'
#' @param capsis_folderpath Folderpath of the capsis source directory
#' @param inv_fp Filepath of the samsaralight inventory
#' @param commandfile_fp Filepath of the command file (if not using rcapsis)
#' @param use_rcapsis If running Samsara with RCapsis package
#' @param java_server_started Boolean, if False, no not start and clode java server 
#'  (to run multiple simulations one after the other without restarting java server)
#'
#' @noRd
sl_run_samsara <- function(capsis_folderpath, inv_fp,
                           use_rcapsis = TRUE,
                           java_server_started = FALSE,
                           commandfile_fp = NULL) {
  
  if (!use_rcapsis) {
    
    if (is.null(commandfile_fp)) {
      stop("No command file provided")
    }
    
    # Capsis command prep
    cmd <- paste("cd", capsis_folderpath,
                 "& capsis -p script samsara2.pgms.flexiblescript.ScriptFlexibleSimuFormat2",
                 commandfile_fp,
                 sep = " ")
    
    # Run command
    tryCatch(shell(cmd), error = function(e) e)
    
    # Create output file
    output <- list("tree" = NULL, "cells" = NULL)
    
    return(output)
  }
  
  # Connect to capsis
  if (!java_server_started) {
    setCapsisPath(path = capsis_folderpath)
    connectToCapsis()
  }
  
  message("Doing analysis...")
  
  # Creates and initializes the script.
  script <- createJavaObject("capsis.app.C4Script", "samsara2")
  
  ip <- createJavaObject("samsara2.model.Samsa2InitialParameters")
  ip$fileName <- inv_fp
  
  script$init(ip)
  
  # Get model and trees
  root <- script$getRoot()
  scene <- root$getScene()
  trees_hashset <- scene$getTrees()
  
  trees <- createJavaObject("java.util.ArrayList")
  trees$addAll(trees_hashset)
  n_trees <- trees$size()
  
  # Get trees energy
  out_trees <- vector(mode = "list", length = n_trees)
  for (i in 0:(n_trees-1)) {
    
    tree <- trees$get(as.integer(i))
    
    out_trees[[i+1]] <- list(
      id_tree = tree$getId(),
      epot = tree$getPotCrownEnergy(),
      e = tree$getEnergy()
    )
  }
  out_trees <- data.table::rbindlist(out_trees)
  
  
  # Get cells energy
  plot <- scene$getPlot()
  cells_hashset <- plot$getCells()
  cells <- createJavaObject("java.util.ArrayList")
  cells$addAll(cells_hashset)
  n_cells<- cells$size()
  
  out_cells <- vector(mode = "list", length = n_cells)
  for (c in 0:(n_cells-1)) {
    
    cell <- cells$get(as.integer(c))
    
    out_cells[[c+1]] <- list(
      x_center = cell$getXCenter(),
      y_center = cell$getYCenter(),
      z_center = cell$getZCenter(),
      e = cell$getTotalEnergy(),
      erel = cell$getRelativeSlopeEnergy()
    )
  }
  out_cells <- data.table::rbindlist(out_cells)
  out_cells[, erel := erel / 100]
  out_cells[order(-y_center, x_center)]
  
  
  # Unconnect from Capsis
  if (!java_server_started) {
    shutdownClient()
  }
  
  # Return list of outputs
  list("trees" = out_trees, "cells" = out_cells)
}
