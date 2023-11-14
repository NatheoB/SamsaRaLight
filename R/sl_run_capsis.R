#' Run SamsaraLight from Capsis
#'
sl_run_capsis <- function(capsis_folderpath, inv_fp,
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
