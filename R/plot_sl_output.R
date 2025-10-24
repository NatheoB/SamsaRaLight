#'  @param plot.trees Character string indicating the filling variables of the trees. 
#'  Set it to NULL if you do not want to plot the trees
#'
#' @import ggplot2 ggforce
#'
#' @export
#'
plot_sl_output <- function(sl_output, 
                           trees.border.species = FALSE,
                           trees.fill = "species",
                           trees.fill.inverse = FALSE,
                           trees.only_inv = FALSE,
                           cells.border = FALSE,
                           cells.fill = NULL,
                           cells.fill.palette = c("base", "base01",
                                                  "viridis", "viridis01",
                                                  "light", "light01"),
                           sensors.plot = FALSE) {
  
  # Check arguments
  if (!cells.fill %in% names(sl_output$output$cells) && !is.null(cells.fill)) {
    stop("cells.fill argument must be a variable within the SamsaRaLight cells output, or be set to NULL...")
  }
  
  cells.fill.palette <- match.arg(cells.fill.palette, 
                                  choices = c("base", "base01",
                                              "viridis", "viridis01",
                                              "light", "light01"), 
                                  several.ok = FALSE)
  
  if (!trees.fill %in% names(sl_output$output$trees) && !trees.fill %in% names(sl_output$input$trees) && !is.null(trees.fill)) {
    stop("trees.fill argument must be a variable within the SamsaRaLight trees input/output, or be set to NULL...")
  }
  
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
    if (cells.fill.palette == "light") {
      plt <- plt +
        scale_fill_gradient(low = grey(0.2), high = grey(1),
                            guide = guide_colorbar(direction = "horizontal"))
    } else if (cells.fill.palette == "light01") {
      plt <- plt +
        scale_fill_gradient(low = grey(0.2), high = grey(1), 
                            limits = c(0, 1),
                            guide = guide_colorbar(direction = "horizontal"))
    } else if (cells.fill.palette == "viridis") {
      plt <- plt +
        scale_fill_viridis_c(guide = guide_colorbar(direction = "horizontal"))
    } else if (cells.fill.palette == "viridis01") {
      plt <- plt +
        scale_fill_viridis_c(guide = guide_colorbar(direction = "horizontal"),
                             limits = c(0, 1))
    } else if (cells.fill.palette == "base") {
      plt <- plt +
        scale_fill_continuous(guide = guide_colorbar(direction = "horizontal"))
    } else if (cells.fill.palette == "base01") {
      plt <- plt +
        scale_fill_continuous(guide = guide_colorbar(direction = "horizontal"),
                              limits = c(0, 1))
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
    if ("added_to_fill" %in% colnames(data_trees) & trees.only_inv) {
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
    } else {
      aes_map_trees <- aes_base_trees
    }
    
    ## Plot trees
    plt <- plt + 
      ggnewscale::new_scale_fill() +
      geom_ellipse(data = data_trees, mapping = aes_map_trees,
                   linewidth = 0.8)
    
    ## Set viridis palette for continuous variables
    if (typeof(data_trees[[trees.fill]]) == "double") {
      
      direction_viridis <- if_else(trees.fill.inverse, -1, 1)
      
      plt <- plt + 
        scale_fill_viridis_c(guide = guide_colorbar(direction = "horizontal"),
                             direction = direction_viridis)
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
