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
#' The plot function also generates a compass indicating:
#' \itemize{
#'   \item Plot orientation (north2x, north in red)
#'   \item Terrain aspect (downslope direction)
#'   \item Slope degrees (annotation only)
#' }
#'
#' @return A `ggplot` object representing the stand.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_curve geom_tile geom_polygon geom_rect labs coord_equal theme_bw theme scale_x_continuous scale_y_continuous xlab ylab facet_wrap annotate scale_colour_manual theme_void geom_point guides
#' @importFrom patchwork plot_layout 
#' @importFrom grid arrow unit
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
      "x" = "X-axis",
      "y" = "Y-axis"
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
          view == "x" ~ x,
          view == "y" ~ y
        ),
        
        # define left/right radius given the view
        r_left_m = case_when(
          view == "x" ~ rxmin_m,
          view == "y" ~ rymin_m
        ),
        r_right_m = case_when(
          view == "x" ~ rxmax_m,
          view == "y" ~ rymax_m
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
      
      labs(x = "Axis position (in m)", y = "Height (m)",
           title = "SamsaraLight input stand",
           subtitle = "top-down view") +
      coord_equal() +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            plot.subtitle = element_text(hjust = 0.5)) +
      facet_wrap(~view_label, ncol = 1)
    
    # Add sensors
    if (add_sensors && !is.null(sl_stand$sensors) && nrow(sl_stand$sensors) > 0) {
      
      sensors_topdown <- sl_stand$sensors %>%
        tidyr::crossing(view = names(view_description)) %>%
        dplyr::mutate(
          pos = case_when(
            view == "x" ~ x,
            view == "y"  ~ y
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
    
    plt_stand <- ggplot() +
      coord_equal() +
      
      # CELLS 
      geom_tile(data = sl_stand$cells, 
                mapping = aes(x = x_center, y = y_center), 
                fill = "white", color = "darkgray")
    
    
    # CORE POLYGON
    if (!is.null(sl_stand$core_polygon$df)) {  
      plt_stand <- plt_stand +
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
      
      plt_stand <- plt_stand +
        geom_ellipse(
          data = trees_plot[i,],
          mapping = aes(
            x0 = x,
            y0 = y,
            a = (rxmax_m + rxmin_m) / 2, # In X-axis
            b = (rymax_m + rymin_m) / 2, # In Y-axis
            angle = 0,
            fill = species
          ),
          color = "black",
          alpha = alpha_val
        )
      
    }
    
    
    # SENSORS
    if (add_sensors && !is.null(sl_stand$sensors) && nrow(sl_stand$sensors) > 0) {
      plt_stand <- plt_stand +
        geom_rect(data = sl_stand$sensors,
                  mapping = aes(xmin = x - 0.5,
                                ymin = y - 0.5,
                                xmax = x + 0.5,
                                ymax = y + 0.5),
                  color = "red", fill = "black")
    }
    
    
    # Axes, titles and theme
    xbreaks <- scales::pretty_breaks(n = 7)(seq(0, sl_stand$geometry$n_cells_x * sl_stand$geometry$cell_size, by = sl_stand$geometry$cell_size))
    ybreaks <- scales::pretty_breaks(n = 7)(seq(0, sl_stand$geometry$n_cells_y * sl_stand$geometry$cell_size, by = sl_stand$geometry$cell_size))
    
    plt_stand <- plt_stand +
      scale_x_continuous(breaks = xbreaks) +
      scale_y_continuous(breaks = ybreaks) +
      xlab("") + ylab("") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            legend.position = "top")
    
    
    # Mini grobs for north2x, aspect and slope
    plt_compass <- plot_orientation_compass(sl_stand$geometry$north2x, 
                                            sl_stand$geometry$slope,
                                            sl_stand$geometry$aspect)
    
    # Create final plot
    plt <- (plt_stand | plt_compass ) +
      plot_layout(widths = c(3, 1),
                  guides = "collect") &
      patchwork::plot_annotation(
        title = "SamsaRaLight input stand",
        subtitle = paste0("\nInventory zone (yellow): ",
                          round(sl_stand$transform$core_area_ha, 2), "ha - ",
                          round(sl_stand$transform$core_batot_m2ha, 2), "m2/ha - ",
                          sum(!x$trees$added_to_fill), " trees",
                          "\nVirtual plot (rectangle):",
                          round(sl_stand$transform$new_area_ha, 2), "ha - ",
                          round(sl_stand$transform$new_batot_m2ha, 2), "m2/ha - ",
                          nrow(x$trees), " trees")
      ) &
      theme(
        legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, face="italic")
      ) &
      guides(colour = "none")
    
  }
  
  plt
}



plot_orientation_compass <- function(north2x, slope, aspect, size = 1) {
  
  ### COMPASS ###
  
  # Initialise cardinal vectors in XY plan when north2x = 0 (i.e. north pointing owards X-axis)
  dirs_north2x_0 <- data.frame(
    dir = c("W", "N", "E", "S"),
    x = c(0, 1, 0, -1),
    y = c(1, 0, -1, 0)
  )
  
  # Rotate compass by north2x
  # north2x is CW from north vector to X axis
  # so if we rotate the X-axis, we rotate CW
  # if we rotate the compass, it is different, we rotate CCW
  theta_north2x <- deg2rad(north2x)
  rot_dirs <- rotate_vec_ccw(dirs_north2x_0$x, dirs_north2x_0$y, theta_north2x)
  
  df_dirs <- data.frame(
    dir = dirs_north2x_0$dir,
    x0 = 0, y0 = 0,
    x1 = size * rot_dirs$x,
    y1 = size * rot_dirs$y,
    col = dirs_north2x_0$dir == "N"
  )
  
  ### ASPECT ###
  
  # Here, we initialise the aspect arrow coordinates to the rotated north compass from above
  df_aspect <- df_dirs[df_dirs$dir == "N",]
  
  # And we rotate by -aspect (because aspect is CW from North)
  theta_aspect <- deg2rad(aspect)
  rot_aspect <- rotate_vec_ccw(df_aspect$x1, df_aspect$y1, -theta_aspect)
  
  df_aspect$x1 <- rot_aspect$x
  df_aspect$y1 <- rot_aspect$y
  
  
  
  ### FINAL PLOT ####
  
  gg <- ggplot() +
    
    # Cardinal directions
    geom_segment(
      data = df_dirs,
      aes(x = x0, y = y0, xend = x1, yend = y1, colour = col),
      arrow = arrow(length = unit(2, "mm")),
      linewidth = 1
    ) +
    
    geom_text(
      data = df_dirs,
      aes(x = 1.2 * x1, y = 1.2 * y1, label = dir),
      size = 4
    ) +
    
    scale_colour_manual(values = c("grey40", "red")) +
    theme_void() +
    theme(legend.position = "none")
  
  # Add aspect arrow ONLY if slope > 0
  if (!is.na(slope) && slope > 0) {
    
    gg <- gg +
      geom_segment(
        data = df_aspect,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        arrow = arrow(length = unit(2, "mm")),
        linewidth = 1,
        colour = "steelblue",
        linetype = "dashed"
      )
  }
  
  gg +
  annotate(
    "text",
    x = 0, y = 1.75 * size,
    label = paste0(
      "slope: ", round(slope, 1), "°"
    ),
    size = 3,
    hjust = 0.5,
    colour = "grey30"
  ) +
    
    coord_equal(
      xlim = c(-1.4, 1.4),
      ylim = c(-1.6, 1.9)
    )
}
