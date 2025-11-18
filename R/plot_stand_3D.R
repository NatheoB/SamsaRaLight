#' Function to plot SamsaRaLight stand in 3D
#' 
#' @param trees a tree table
#'
#' @import plotly
#'
#' @export
#'
plot_stand_3D <- function(trees = NULL,
                          species_colors = NULL,
                          slope = 0, 
                          aspect = 0, 
                          north_to_x_cw = 90) {
  
  if (any(!unique(trees$species) %in% names(species_colors))) {
    stop("Some species do not have specified colors")
  }
  
  # FUNCTIONS
  make_mesh_E <- function(x, y, z_ground, hbase, htop, radius,
                          n_theta = 30, n_phi = 15) {
    
    # Compute the ellispoid parameters
    a <- radius
    b <- radius
    
    # Vertical radius and center
    c <- (htop - hbase) / 2
    z0 <- (htop + hbase) / 2
    
    # Create theta and phi grid
    theta <- seq(0, pi, length.out = n_theta)
    phi <- seq(0, 2*pi, length.out = n_phi)
    theta_grid <- matrix(rep(theta, each=length(phi)), nrow=length(phi))
    phi_grid <- matrix(rep(phi, length(theta)), nrow=length(phi))
    
    # Ellipsoid coordinates
    x <- a * sin(theta_grid) * cos(phi_grid) + x
    y <- b * sin(theta_grid) * sin(phi_grid) + y
    z <- c * cos(theta_grid) + z0 + z_ground
    
    list(x = x, y = y, z = z) 
  }
  
  make_mesh_P <- function(x, y, z_ground, hbase, htop, radius,
                          n_theta = 30, n_phi = 15) {
    
    # Create theta and phi grid
    theta <- seq(0, pi/2, length.out = n_theta)  # only 0 to pi/2 for paraboloid (top to bottom)
    phi <- seq(0, 2*pi, length.out = n_phi)
    
    theta_grid <- matrix(rep(theta, each=length(phi)), nrow=length(phi))
    phi_grid <- matrix(rep(phi, length(theta)), nrow=length(phi))
    
    # Radial distance from center
    r <- radius * sin(theta_grid)  # sin(theta) to go from 0 (top) to radius (base)
    
    # Paraboloid height
    z <- htop - (htop - hbase) * (r^2 / radius^2)
    
    # Convert to Cartesian coordinates
    x <- r * cos(phi_grid) + x
    y <- r * sin(phi_grid) + y
    z <- z + z_ground
    
    list(x = x, y = y, z = z)
  }
  
  
  # Compute bottom azimut of the stand
  bottom_azimut <- get_bottom_azimut(aspect, north2x)
  
  # Compute tree z coordinate
  trees$z_ground <- -get_z(trees$x, trees$y,
                           deg2rad(slope), deg2rad(bottom_azimut))
  
  
  # Init the plot
  plt <- plot_ly()
  
  
  # ---- Add 3D crowns (ellipsoids) ----
  for(i in 1:nrow(trees)) {
    
    if (trees$crown_type[i] %in% c("E", "2E", "8E")) {
      crown_mesh <- make_mesh_E(
        x = trees$x[i],
        y = trees$y[i],
        z_ground = trees$z_ground[i],
        hbase = trees$hbase_m[i],
        htop  = trees$h_m[i],
        radius = mean(c(trees$rn_m[i], trees$rs_m[i], trees$re_m[i], trees$rw_m[i])),
        n_theta = 5, 
        n_phi = 5
      )
    } else if (trees$crown_type[i] %in% c("P", "4P")) {
      crown_mesh <- make_mesh_P(
        x = trees$x[i],
        y = trees$y[i],
        z_ground = trees$z_ground[i],
        hbase = trees$hbase_m[i],
        htop  = trees$h_m[i],
        radius = mean(c(trees$rn_m[i], trees$rs_m[i], trees$re_m[i], trees$rw_m[i])),
        n_theta = 5,
        n_phi = 5
      )
    } else {
      stop("Wrong crown type")
    }
    
    # Unique surfacecolor per tree
    surface_color_matrix <- matrix(1, nrow=nrow(crown_mesh$x), ncol=ncol(crown_mesh$x))
    
    # Add the surface
    plt <- plt %>%
      add_surface(
        x = crown_mesh$x,
        y = crown_mesh$y,
        z = crown_mesh$z,
        surfacecolor = surface_color_matrix,
        colorscale = list(c(0,1), c(species_colors[trees$species[i]], 
                                    species_colors[trees$species[i]])),
        cmin = 0, 
        cmax = 1,
        showscale = FALSE,
        opacity = 0.7
      )
  }
  
  
  # ---- Add trunks ----
  df_trunks <- do.call(rbind, lapply(1:nrow(trees), function(i) {
    data.frame(
      x = c(trees$x[i], trees$x[i]),
      y = c(trees$y[i], trees$y[i]),
      z = c(trees$z_ground[i], trees$z_ground[i] + trees$h_m[i]),
      id = trees$id_tree[i]
    )
  }))
  
  plt <- plt %>%
    add_trace(
      data = df_trunks,
      x = ~x, y = ~y, z = ~z,
      split = ~id,
      type = "scatter3d",
      mode = "lines",
      line = list(width = 3, color = "sienna"), 
      inherit = FALSE
    )
  
  
  # ---- Terrain ----
  df_geom_terrain <- data.frame(
    x = rep(c(0, size_x), each = 2),
    y = rep(c(0, size_y), time = 2)) %>%
    dplyr::mutate(z = -get_z(x, y,
                             deg2rad(slope),
                             deg2rad(bottom_azimut)))
  
  # ---- Center of the terrain ----
  size_x = max(trees$x) - min(trees$x)
  size_y = max(trees$y) - min(trees$y)
  
  cx <- (max(trees$x) - min(trees$x)) / 2
  cy <- (max(trees$y) - min(trees$y)) / 2
  
  
  # ---- Cardinal directions ----
  theta <- deg2rad(aspect - north2x)
  r <- min(size_x, size_y) / 2
  
  
  # East–West line
  df_geom_lineEW <- data.frame(
    x = cx + c( r * cos(theta + pi/2),
                r * cos(theta + 3*pi/2)),
    y = cy + c(-r * sin(theta + pi/2),
               -r * sin(theta + 3*pi/2)),
    z = 0,
    lab = c("E", "W")
  )
  
  # North–South line
  df_geom_lineNS <- data.frame(
    x = cx + c( r * cos(theta),
                r * cos(theta + pi)),
    y = cy + c(-r * sin(theta),
               -r * sin(theta + pi)),
    z = 0,
    lab = c("N", "S")
  )
  
  
  plt <- plt %>%
    
    # Plot the compass
    add_trace(data = df_geom_lineNS,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(width = 2, color = "red")) %>%
    add_text(data = df_geom_lineNS,
             x = ~x, y = ~y, z = ~z,
             text = ~lab, mode = "markers",
             textfont = list(size = 12, color = "black")) %>%
    
    add_trace(data = df_geom_lineEW,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(width = 2, color = "red")) %>%
    add_text(data = df_geom_lineEW,
             x = ~x, y = ~y, z = ~z,
             text = ~lab, mode = "markers",
             textfont = list(size = 12, color = "black")) %>%
    
    # Plot the terrain
    add_mesh(data = df_geom_terrain,
             x = ~x, y = ~y, z = ~z,
             intensity = ~z,
             opacity = 1,
             colorscale = 'Viridis') %>%
    
    hide_guides() %>% 
    layout(
      scene = list(
        # camera position and target
        camera = list(
          eye = list(x = -1.8, y = -1.8, z = 1.5),
          center = list(x = 0.15, y = 0.15, z = -0.1),
          up = list(x = 0, y = 0, z = 1)
        ),
        
        # automatically set axis ranges (to fit terrain)
        xaxis = list(autorange = TRUE),
        yaxis = list(autorange = TRUE),
        zaxis = list(autorange = TRUE)
      )
    )
  
  print(plt)

}
