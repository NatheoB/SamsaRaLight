#' Function to plot SamsaRaLight object output
#' 
#'  @param plot.trees Character string indicating the filling variables of the trees. 
#'  Set it to NULL if you do not want to plot the trees
#'
#' @import ggplot2 ggforce
#'
#' @export
#'
plot_sl_output <- function(sl_output, 
                           trees.only_inventoried = FALSE,
                           trees.border.species = FALSE,
                           trees.fill = "e",
                           trees.fill.palette = c("viridis", "light", "base"),
                           trees.fill.inverse = FALSE,
                           trees.fill.limits = NULL,
                           cells.border = FALSE,
                           cells.fill = "e",
                           cells.fill.palette = c("light", "viridis", "base"),
                           cells.fill.inverse = FALSE,
                           cells.fill.limits = NULL,
                           sensors.plot = FALSE) {
  
  # Check arguments
  if (!cells.fill %in% names(sl_output$output$cells) && !is.null(cells.fill)) {
    stop("cells.fill argument must be a variable within the SamsaRaLight cells output, or be set to NULL...")
  }
  
  cells.fill.palette <- match.arg(cells.fill.palette, 
                                  choices = c("light", "viridis", "base"), 
                                  several.ok = FALSE)
  
  if (!trees.fill %in% names(sl_output$output$trees) && !trees.fill %in% names(sl_output$input$trees) && !is.null(trees.fill)) {
    stop("trees.fill argument must be a variable within the SamsaRaLight trees input/output, or be set to NULL...")
  }
  
  trees.fill.palette <- match.arg(trees.fill.palette, 
                                  choices = c("viridis", "light", "base"), 
                                  several.ok = FALSE)
  
  # Prepare the plot grid
  cell_size <- sl_output$input$info$cell_size
  n_cells_x <- sl_output$input$info$n_cells_x
  n_cells_y <- sl_output$input$info$n_cells_y
  
  # Base plot
  plt <- ggplot() +

    coord_equal() +
    scale_x_continuous(breaks = seq(0, n_cells_x*cell_size,
                                    by = cell_size),
                       labels = round(seq(0, n_cells_x*cell_size,
                                          by = cell_size),
                                      digits = 1)) +
    scale_y_continuous(breaks = seq(0, n_cells_y*cell_size,
                                    by = cell_size),
                       labels = round(seq(0, n_cells_y*cell_size,
                                          by = cell_size),
                                      digits = 1)) +
    
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right") +
    xlab("") + ylab("")
  
  
  # Fill cells if precised  
  if (!is.null(cells.fill)) {

    plt <- plt + 
      ggnewscale::new_scale_fill() +
      geom_tile(data = sl_output$output$cells, 
                mapping = aes(
                  x = x_center,
                  y = y_center,
                  fill = !!sym(cells.fill)
                ),
                color = if_else(cells.border, "black", NA))
    
    
    ## Set palette for cells
    guide_cells <- guide_colorbar(direction = "horizontal",
                                  title.position = "top",
                                  title.hjust = 0.5)
    
    if (cells.fill.palette == "light") {
      
      col_low_cells <- if_else(cells.fill.inverse, grey(1), grey(0.2))
      col_high_cells <- if_else(cells.fill.inverse, grey(0.2), grey(1))
      
      
      plt <- plt +
        scale_fill_gradient(low = col_low_cells, high = col_high_cells,
                            guide = guide_cells,
                            limits = cells.fill.limits)
      
    } else if (cells.fill.palette == "viridis") {
      
      direction_color_cells <- if_else(cells.fill.inverse, -1, 1)
      
      plt <- plt +
        scale_fill_viridis_c(guide = guide_cells,
                             direction = direction_color_cells,
                             limits = cells.fill.limits)
      
    } else if (cells.fill.palette == "base") {
      
      col_low_cells <- if_else(cells.fill.inverse, "#56B1F7", "#132B43")
      col_high_cells <- if_else(cells.fill.inverse, "#132B43", "#56B1F7")
      
      plt <- plt +
        scale_fill_continuous(low = col_low_cells, high = col_high_cells,
                              guide = guide_cells,
                              limits = cells.fill.limits)
      
    } 
    

  } else if (cells.border) {
    
    plt <- plt + 
      geom_tile(data = sl_output$output$cells, 
                mapping = aes(
                  x = x_center,
                  y = y_center
                ),
                color = "black", fill = NA)
    
  }
  
  # Add plot border
  plt <- plt +
    geom_rect(data = data.frame(xmin = 0, xmax = n_cells_x*cell_size,
                                ymin = 0, ymax = n_cells_y*cell_size),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              color = "black", fill = NA,
              linewidth = 1.05)
  
  
  
  # Plot trees if precised  
  if (!is.null(trees.fill)) {
    
    data_trees <- dplyr::left_join(
      sl_output$input$trees,
      sl_output$output$trees,
      by = c("id_tree", "x", "y")
    )
    
    ## Remove trees added to fill around the inventory zone if precsied
    ## Only if the column exists, i.e. only if the virtual stand has been created with the create_rect_stand() function
    if ("added_to_fill" %in% colnames(data_trees) & trees.only_inventoried) {
      data_trees <- data_trees %>% 
        dplyr::filter(!added_to_fill)
    }
    
    ## Base mapping (always used)
    aes_base_trees <- aes(
      x0 = x,
      y0 = y,
      a = (re_m + rw_m) / 2,
      b = (rn_m + rs_m) / 2,
      angle = 0,
      fill = !!sym(trees.fill)
    )
    
    ## Add color mapping conditionally
    if (trees.border.species) {
      aes_map_trees <- modifyList(aes_base_trees, aes(color = species))
      trees_linewidth <- 0.8
    } else {
      aes_map_trees <- aes_base_trees
      trees_linewidth <- 0.1
    }
    
    ## Plot trees
    plt <- plt + 
      ggnewscale::new_scale_fill() +
      geom_ellipse(data = data_trees, mapping = aes_map_trees,
                   linewidth = trees_linewidth)
    
    ## Set colors for continuous variables
    if (typeof(data_trees[[trees.fill]]) %in% c("integer", "numeric", "double")) {
      
      guide_trees <- guide_colorbar(direction = "horizontal",
                                    title.position = "top",
                                    title.hjust = 0.5)
      
      if (trees.fill.palette == "light") {
        
        col_low_trees <- if_else(trees.fill.inverse, grey(1), grey(0.2))
        col_high_trees <- if_else(trees.fill.inverse, grey(0.2), grey(1))
        
        
        plt <- plt +
          scale_fill_gradient(low = col_low_trees, high = col_high_trees,
                              guide = guide_trees,
                              limits = trees.fill.limits)
        
      } else if (trees.fill.palette == "viridis") {
        
        direction_color_trees <- if_else(trees.fill.inverse, -1, 1)
        
        plt <- plt +
          scale_fill_viridis_c(guide = guide_trees,
                               direction = direction_color_trees,
                               limits = trees.fill.limits)
        
      } else if (trees.fill.palette == "base") {
        
        col_low_trees <- if_else(trees.fill.inverse, "#56B1F7", "#132B43")
        col_high_trees <- if_else(trees.fill.inverse, "#132B43", "#56B1F7")
        
        plt <- plt +
          scale_fill_continuous(low = col_low_trees, high = col_high_trees,
                                guide = guide_trees,
                                limits = trees.fill.limits)
        
      } 
      
    }
  
  }
  
  # Plot sensors if specified
  if (sensors.plot & nrow(sl_output$output$sensors) > 0) {

    plt <- plt +
      geom_rect(data = sl_output$output$sensors,
                mapping = aes(xmin = x - 1,
                              ymin = y - 1,
                              xmax = x + 1,
                              ymax = y + 1),
                color = "red", fill = "black")
    
  }
  
  # Return the plot
  plt
}
