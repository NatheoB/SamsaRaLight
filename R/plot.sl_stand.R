#' Plot a SamsaRaLight virtual stand
#'
#' This function plots a virtual forest stand (`sl_stand`) produced by SamsaRaLight.
#' It can display a top-down view with tree crowns or a side/top view with cells and trees.
#'
#' @param x An object of class `sl_stand`.
#' @param ... Additional arguments passed to lower-level plotting functions.
#' @param top_down Logical, if TRUE, creates a top-down view with multiple directions (south, north, west, east).
#' @param only_inv Logical, if TRUE, plot only trees from the initial inventory (i.e. not trees added to fill around the core polygon)
#' @param add_sensors Logical; if TRUE (default), sensors are drawn on the plot.
#'   In top-down mode, sensors are shown as segment from ground to their height; in map view,
#'   sensors are drawn as squares on the ground.
#'
#' @details
#' For the sake of the representation in top-down plot, z are offset such as minimum altitude tree is at Y-axis height = 0
#' 
#'
#' @return A `ggplot` object representing the stand.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_curve geom_tile geom_polygon geom_rect labs coord_equal theme_bw theme scale_x_continuous scale_y_continuous xlab ylab facet_wrap guides guide_legend
#' @importFrom ggforce geom_ellipse
#' @importFrom dplyr filter mutate case_when
#' @importFrom tidyr crossing
#' 
#' @export
#' @method plot sl_stand
#' 
plot.sl_stand <- function(x, ..., 
                          top_down = FALSE,
                          only_inv = FALSE,
                          add_sensors = TRUE) {
  
  stopifnot(inherits(x, "sl_stand"))
  stopifnot(is.logical(top_down), length(top_down) == 1)
  stopifnot(is.logical(only_inv), length(only_inv) == 1)
  
  sl_stand <- x  # Rename for internal use
  
  if (only_inv) {
    sl_stand$trees <- sl_stand$trees %>% 
      dplyr::filter(!added_to_fill)
  }
  
  plt <- NULL
  
  ## Top-down view
  if (top_down) {
    
    # Label for plots 
    view_description <- c(
      "south" = "south view (X-axis W->E)",
      # "west" = "west view (Y-axis N->S)",
      # "north" = "north view (X-axis E->W)",
      "east" = "east view (Y-axis S->N)"
    )
    
    # z-offset
    z_offset <- - min(sl_stand$trees$z)
    
    # Create the trees for each top-down view
    trees_topdown <- sl_stand$trees %>% 
      tidyr::crossing(view = names(view_description)) %>% 
      dplyr::mutate(
        
        view_label = unname(view_description[view]),
        view_label = factor(view_label, levels = unname(view_description)),
        
        # Change z for better plots (positive altitude)
        # Lowest tree altitude at 0
        z = z + z_offset,
        
        # Height of the maximum radius
        hmax_m = case_when(
          crown_type == "E" ~ (h_m + hbase_m) / 2, 
          crown_type %in% c("P", "4P") ~ hbase_m,
          crown_type %in% c("2E", "8E") ~ hmax_m
        ),
        
        # define tree position given the view
        pos = case_when(
          view == "south" ~ x,
          view == "north" ~ - x,
          view == "west" ~ - y,
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
      
      labs(x = "Axis position (in m)", y = "Height (m)") +
      coord_equal() +
      theme_bw() +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
      facet_wrap(~view_label)
    
    # Add sensors
    if (add_sensors && !is.null(sl_stand$sensors) && nrow(sl_stand$sensors) > 0) {
      
      sensors_topdown <- sl_stand$sensors %>%
        tidyr::crossing(view = names(view_description)) %>%
        dplyr::mutate(
          pos = case_when(
            view == "south" ~ x,
            # view == "north" ~ - x,
            # view == "west"  ~ - y,
            view == "east"  ~ y
          ),
          z = z + z_offset,
          view_label = unname(view_description[view]),
          view_label = factor(view_label, levels = unname(view_description))
        )
      
      plt <- plt +
        geom_segment(data = sensors_topdown,
                     aes(x = pos, xend = pos, 
                         y = z - h_m, yend = z), 
                     linewidth = 1,
                     color = "red")
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
    
    # Set alpha for transparency
    alpha_val <- 0.6
    
    # Order trees by height (smaller first, taller on top)
    trees_plot <- sl_stand$trees[order(sl_stand$trees$h_m), ]
    
    for (i in seq_len(nrow(trees_plot))) {
      
      plt <- plt +
        geom_ellipse(
          data = trees_plot[i,],
          mapping = aes(
            x0 = x,
            y0 = y,
            a = (re_m + rw_m) / 2,
            b = (rn_m + rs_m) / 2,
            angle = 0,
            fill = species
          ),
          color = "black",
          alpha = alpha_val
        )
      
    }
    
    
    # SENSORS
    if (add_sensors && !is.null(sl_stand$sensors) && nrow(sl_stand$sensors) > 0) {
      plt <- plt +
        geom_rect(data = sl_stand$sensors,
                  mapping = aes(xmin = x - 0.5,
                                ymin = y - 0.5,
                                xmax = x + 0.5,
                                ymax = y + 0.5),
                  color = "red", fill = "black")
    }
    
    
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
      
      labs(title = "SamsaRaLight input stand",
           subtitle = paste0("\nInventory zone (yellow): ",
                             round(sl_stand$transform$core_area_ha, 2), "ha - ",
                             round(sl_stand$transform$core_batot_m2ha, 2), "m2/ha - ",
                             sum(!x$trees$added_to_fill), " trees",
                             "\n\nVirtual plot (rectangle): ",
                             round(sl_stand$transform$new_area_ha, 2), "ha - ",
                             round(sl_stand$transform$new_batot_m2ha, 2), "m2/ha - ",
                             nrow(x$trees), " trees")) +
      
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) 
    
    
  }
  
  plt
}