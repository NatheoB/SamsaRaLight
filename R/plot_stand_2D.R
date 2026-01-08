#' Function to plot SamsaRaLight trees with top/down view
#' 
#' @param trees a tree table
#'
#' @import ggplot2 ggforce cowplot
#'
#' @export
#'
plot_stand_2D <- function(sl_stand, 
                          top_down = FALSE, 
                          transparency = TRUE) {
  
  
  plt <- NULL
  
  ## Top-down view
  if (top_down) {
    
    # Label for plots 
    view_description <- c(
      "south" = "south view - W to E",
      "north" = "north view - E to W",
      "west" = "west view - N to S",
      "east" = "east view - S to N"
    )
    
    
    # Create the plot for each top-down view
    for (view in names(view_description)) {
      
      trees_topdown <- sl_stand$trees %>% 
        dplyr::mutate(
          
          view_label = view_description[[view]],
          
          # define tree position given the view
          pos = case_when(
            view == "south" ~ x,
            view == "north" ~ -x,
            view == "west" ~ -y,
            view == "east" ~ y
          ),
          
          # define left/right radius given the view
          r_left_m = case_when(
            view == "south" ~ rw_m,
            view == "north" ~ re_m,
            view == "west" ~ rn_m,
            view == "east" ~ rs_m
          ),
          r_right_m = case_when(
            view == "south" ~ re_m,
            view == "north" ~ rw_m,
            view == "west" ~ rs_m,
            view == "east" ~ rn_m
          )
        )
      
      
      # Order trees 
      if (view == "south") {
        trees_topdown %>% dplyr::arrange(desc(y))
      } else if (view == "north") {
        trees_topdown %>% dplyr::arrange(y)
      } else if (view == "west") {
        trees_topdown %>% dplyr::arrange(desc(x))
      } else if (view == "east") {
        trees_topdown %>% dplyr::arrange(x)
      }
      
      # 2D top/down
      plt <- ggplot(trees_topdown) +
        
        # trunks
        geom_segment(aes(x = pos, xend = pos, 
                         y = z, yend = z + hmax_m), 
                     linewidth = 0.1) +
        
        # --- Ellipsoids ---
        geom_curve(
          data = subset(trees_topdown, crown_type %in% c("E", "2E", "8E")),
          aes(x = pos - r_left_m, xend = pos, 
              y = z + hmax_m, yend = z + h_m, 
              color = species),
          curvature = -0.4
        ) +
        geom_curve(
          data = subset(trees_topdown, crown_type %in% c("E", "2E", "8E")),
          aes(x = pos + r_right_m, xend = pos, 
              y = z + hmax_m, yend = z + h_m, 
              color = species),
          curvature = 0.4
        ) +
        
        # --- Paraboloids ---
        geom_curve(
          data = subset(trees_topdown, crown_type %in% c("P", "4P")),
          aes(x = pos - r_left_m, xend = pos, 
              y = z + hmax_m, yend = z + h_m, 
              color = species),
          curvature = -0.2
        ) +
        geom_curve(
          data = subset(trees_topdown, crown_type %in% c("P", "4P")),
          aes(x = pos + r_right_m, xend = pos, 
              y = z + hmax_m, yend = z + h_m, 
              color = species),
          curvature = 0.2
        ) +
        
        labs(x = "Position", y = "Height (m)") +
        coord_equal() +
        theme_bw() +
        facet_wrap(~view_label, ncol = 1)
      
    }
      
  }
  
  # Plot from above view
  else {
    
    plt <- ggplot() +
      coord_equal() +
      
      # CELLS 
      geom_tile(data = sl_stand$cells, 
                mapping = aes(x = x_center, y = y_center), 
                fill = "white", color = "darkgray")
    
    
    # CORE POLYGON
    if (!is.null(sl_stand$core_polygon$df)) {  
      plt <- plt +
        geom_polygon(data = sl_stand$core_polygon$df,
                     mapping = aes(x = x, y = y),
                     fill = "yellow", color = "black", alpha = 0.5)
    }
    
    
    # TREES
    
    # Set alpha based on transparency argument
    alpha_val <- ifelse(transparency, 0.6, 1)
    
    # Order trees by height (smaller first, taller on top)
    trees_plot <- sl_stand$trees[order(sl_stand$trees$h_m), ]
    
    # Base aesthetics
    aes_mapping_trees <- aes(
      x0 = x,
      y0 = y,
      a = (re_m + rw_m) / 2,
      b = (rn_m + rs_m) / 2,
      angle = 0
    )
    
    # Add species fill if available
    if ("species" %in% names(trees_plot)) {
      aes_mapping_trees$fill <- quote(species)
    }
    
    for (i in seq_len(nrow(trees_plot))) {
      
      if ("species" %in% names(trees_plot)) {
        plt <- plt +
          geom_ellipse(
            data = trees_plot[i,],
            mapping = aes_mapping_trees,
            color = "black",
            alpha = alpha_val
          )
      } else {
        plt <- plt +
          geom_ellipse(
            data = trees_plot[i,],
            mapping = aes_mapping_trees,
            color = "black",
            alpha = alpha_val,
            fill = "grey"
          )
      }
      
    }
    
    
    # SENSORS
    # plt <- plt +
    #   geom_rect(data = data_stand$sensors,
    #             mapping = aes(xmin = x - 0.5,
    #                           ymin = y - 0.5,
    #                           xmax = x + 0.5,
    #                           ymax = y + 0.5),
    #             color = "red", fill = "black")
      
      # GRAPHIC
    plt <- plt +
      scale_x_continuous(breaks = seq(0, sl_stand$geometry$n_cells_x * sl_stand$geometry$cell_size,
                                      by = sl_stand$geometry$cell_size),
                         labels = round(seq(0, sl_stand$geometry$n_cells_x * sl_stand$geometry$cell_size,
                                            by = sl_stand$geometry$cell_size),
                                        digits = 1)) +
      scale_y_continuous(breaks = seq(0, sl_stand$geometry$n_cells_y * sl_stand$geometry$cell_size,
                                      by = sl_stand$geometry$cell_size),
                         labels = round(seq(0, sl_stand$geometry$n_cells_y * sl_stand$geometry$cell_size,
                                            by = sl_stand$geometry$cell_size),
                                        digits = 1)) +
      xlab("") + ylab("") +
      
      labs(title = "SamsaRaLight virtual stand",
           subtitle = paste0(round(sl_stand$transform$new_area_ha, 2), "ha - ",
                             round(sl_stand$transform$new_batot_m2ha, 2), "m2/ha - ",
                             nrow(sl_stand$trees), " trees")) +
      
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) 
    
    
    # Add terrain info
  }
  
  plt
}