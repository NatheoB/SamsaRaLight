#' Function to plot SamsaRaLight trees with top/down view
#' 
#' @param trees a tree table
#'
#' @import ggplot2 ggforce cowplot
#'
#' @export
#'
plot_stand_2D <- function(trees, 
                          top_down = FALSE,
                          slope = 0, 
                          aspect = 0, 
                          north_to_x_cw = 90) {
  
  
  if (top_down) {
    
    # Get stand geometry variables
    slope_rad <- deg2rad(slope)
    bottom_azimuth_rad <- deg2rad(get_bottom_azimut(aspect, north_to_x_cw))
    
    # Compute z position
    trees <- trees %>% 
      dplyr::mutate(
        z = get_z(x, y, slope_rad, bottom_azimuth_rad)
      )
    
    # Label for plots 
    view_description <- c(
      "south" = "south view - W to E",
      "north" = "north view - E to W",
      "west" = "west view - N to S",
      "east" = "east view - S to N"
    )
    
    # Create the trees for each top-down view
    trees_topdown <- trees %>% 
      tidyr::crossing(view = names(view_description)) %>% 
      dplyr::mutate(
        
        view_label = view_description[view],
        
        # Height of the maximum radius
        hmax_m = case_when(
          crown_type == "E" ~ (h_m + hbase_m) / 2, 
          crown_type == "P" ~ hbase_m,
          crown_type == "8E" ~ hmax_m
        ),
        
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
  
  else {
    
    plt <- ggplot() +
      
      geom_ellipse(data = trees, 
                   mapping = aes(
                     x0 = x,
                     y0 = y,
                     a = (re_m + rw_m) / 2,
                     b = (rn_m + rs_m) / 2,
                     angle = 0,
                     fill = species
                   )) +
      
      coord_equal() +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right") +
      xlab("") + ylab("")
    
    # Add geometry information
    plt <- add_plane_geometry(
      main_plot = plt,
      slope  = slope,      # degrees
      aspect = aspect      # degrees clockwise from North
    )
    
  }
  
  print(plt)
}


#-----------------------------------------
# Small compass with aspect arrow
#-----------------------------------------
compass_grob <- function(aspect_deg = 0, size = 1) {
  rad <- aspect_deg * pi/180
  
  ggplot() +
    # Main aspect arrow
    geom_segment(
      aes(x = 0, y = 0,
          xend = size * sin(rad),
          yend = size * cos(rad)),
      arrow = arrow(length = unit(0.25, "cm")),
      linewidth = 1
    ) +
    # Cardinal directions
    geom_segment(aes(x=0, y=0, xend=0, yend=size)) +
    geom_segment(aes(x=0, y=0, xend=size, yend=0)) +
    geom_segment(aes(x=0, y=0, xend=-size, yend=0)) +
    geom_segment(aes(x=0, y=0, xend=0, yend=-size)) +
    geom_text(aes(0, size*1.2, label="N"), size=3.2) +
    geom_text(aes(size*1.2, 0, label="E"), size=3.2) +
    geom_text(aes(0, -size*1.2, label="S"), size=3.2) +
    geom_text(aes(-size*1.2, 0, label="W"), size=3.2) +
    coord_fixed() +
    theme_void()
}

#-----------------------------------------
# Slope indicator (tilted surface + angle)
#-----------------------------------------
slope_grob <- function(slope_deg = 10, aspect_deg = 0, size = 1){
  
  z_drop <- size * tan(slope_deg * pi/180)
  
  ggplot() +
    # Tilted slope
    geom_segment(
      aes(x = 0, y = z_drop,
          xend = size, yend = 0),
      linewidth = 1.2
    ) +
    # Angle label
    geom_text(
      aes(size*0.1, z_drop*0.9,
          label = paste0(slope_deg, "°")),
      size = 3.3
    ) +
    # Down-slope arrow (direction)
    geom_segment(
      aes(x = size*0.5, y = z_drop*0.5,
          xend = size*0.5 + 0.3 * sin(aspect_deg*pi/180),
          yend = z_drop*0.5 + 0.3 * cos(aspect_deg*pi/180)),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "blue"
    ) +
    coord_fixed() +
    theme_void()
}

#-----------------------------------------
# Combine compass + slope and inset in plot
# using cowplot::ggdraw / draw_plot
#-----------------------------------------
add_plane_geometry <- function(main_plot, slope, aspect,
                               compass_size = 1, slope_size = 1,
                               x = 0.72, y = 0.72, # NPC coords (0-1), bottom-left of inset
                               width = 0.25, height = 0.25,
                               align = "hv") {
  
  # Create the small figure (stacked vertically)
  mini_fig <- cowplot::plot_grid(
    compass_grob(aspect_deg = aspect, size = compass_size),
    slope_grob(slope_deg = slope, aspect_deg = aspect, size = slope_size),
    ncol = 1,
    rel_heights = c(1.2, 1)
  )
  
  # Compose: main plot + inset
  out <- cowplot::ggdraw() +
    cowplot::draw_plot(main_plot, 0, 0, 1, 1) +
    cowplot::draw_plot(mini_fig, x = x, y = y, width = width, height = height, scale = 1)
  
  out
}