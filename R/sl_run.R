#' Compute samsara ligth radiative balance
#'
#' @param trees a data.frame with one row for each tree and 9 columns:
#' \itemize{
#'  \item{"id": }{Unique id of the tree}
#'  \item{"x": }{X position of the tree within the stand}
#'  \item{"y": }{Y position of the tree within the stand}
#'  \item{"z": }{Z position of the tree within the stand}
#'  \item{"height": }{Height of the trees (in meters)}
#'  \item{"cradius": }{Crown radius of the tree (in meters)}
#'  \item{"cdepth": }{Crown depth (i.e. height of the crown) of the tree (in meters)}
#'  \item{"crown_lad"}{Leaf Area Density (in m2 of leaves per m3 of crown)}
#'  \item{"crown_p"}{Crown openness (no unit)}
#'  \item{"crown_type"}{Between "P" for paraboloid and "E" for ellipsoid}
#'  \item{"id_cell"}{Unique id of the cell the tree belong to (i.e. `cells` data.frame parameter)}
#' }
#' @param monthly_rad data.frame - Monthly horizontal radiation (Hrad) and diffuse to global ratio (DGratio)
#'    Computed with function samsaRa::sl_get_monthlyrad()
#' @param latitude double - Latitude of the plot (in degrees)
#' @param start_day integer between 1 and 365 - First day of the vegetative period
#' @param end_day integer between 1 and 365 - Last day of the vegetative period
#' @param slope double - Slope of the plot (in degrees)
#' @param north_to_x_cw double - Angle from North to x axis clockwise. (in degrees)
#'    Default correspond to a Y axis oriented toward the North.
#' @param aspect double - Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
#'    northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
#' @param cell_size Length of the side of a squared cell composing the stand (in meters)
#' @param n_cells integer - Number of cells of the side of a squared stand
#' @param use_torus if True, compute light competition using borders modelled with a
#'  torus system, otherwise,borders are open grasslands
#' @param turbid_medium If TRUE, crown are considered as turbid medium, otherwise,
#'  considered as porous envelope
#' @param use_rcpp If TRUE, use Rcpp package and C++ implementation of SamsaraLight for faster
#'  computation. Otherwise, run SamsaraLight with only R-based scripts
#' @param trunk_interception Consider interception of rays by trunks
#'
#' @import data.table
#' @importFrom Rcpp sourceCpp
#' @useDynLib SamsaRaLight, .registration = TRUE
#'
#' @export
#'
sl_run <- function(trees,
                   monthly_rad,
                   latitude = 46,
                   start_day = 1, end_day = 365,
                   slope = 0,
                   north_to_x_cw = 90,
                   aspect = 0,
                   cell_size = 10,
                   n_cells = 10,
                   use_torus = TRUE,
                   turbid_medium = TRUE,
                   use_rcpp = TRUE,
                   trunk_interception = FALSE) {

  
  # Create monthly rays
  rays <- sl_create_monthly_rays(monthly_rad = monthly_rad,
                                 latitude = latitude,
                                 start_day = start_day,
                                 end_day = end_day,
                                 soc = TRUE,
                                 slope = slope,
                                 north_to_x_cw = north_to_x_cw,
                                 aspect = aspect,
                                 height_anglemin = 10,
                                 direct_startoffset = 0,
                                 direct_anglestep = 5,
                                 diffuse_anglestep = 15)

  # Create cells dataframe
  cells <- data.table(
    x_id = rep(0:(n_cells-1), times = n_cells),
    y_id = rep(0:(n_cells-1), each = n_cells)
  )
  
  cells <- cells[, `:=`(x_center = x_id * cell_size + cell_size / 2,
                        y_center = n_cells * cell_size - (y_id * cell_size + cell_size / 2))]
  cells <- cells[, z_center := get_z(x_center, y_center, deg2rad(slope), deg2rad(north_to_x_cw - aspect))]
  
  setorder(cells, x_id, y_id)
  cells <- data.table(id_cell = 1:nrow(cells), cells)
  
  # Search for cell in which the tree belong to
  trees <- as.data.table(trees)
  
  trees <- trees[, `:=`(z = get_z(x, y,
                                  deg2rad(slope),
                                  deg2rad(-aspect + north_to_x_cw)),
                        xid_cell = x %/% cell_size,
                        yid_cell = n_cells - y %/% cell_size - 1)]
  
  trees <- cells[, .(id_cell, x_id, y_id)][trees, on = c(x_id = "xid_cell", y_id = "yid_cell")]
  trees <- trees[, `:=`(x_id = NULL, y_id = NULL)]
  

  # USING RCPP (C++ script in src folder)
  if (use_rcpp) {

    # Run call to c++ script
    out <- sl_run_rcpp(
      trees, cells, rays$rays,
      sum(rays$e_slope),
      slope, north_to_x_cw, aspect,
      cell_size, n_cells,
      use_torus, turbid_medium, trunk_interception)

    # Convert dataframe into data.table
    out$trees <- as.data.table(out$trees)
    out$cells <- as.data.table(out$cells)
    
    return(out)
  }
  
  # OTHERWISE, USE R
  
  slope_rad <- deg2rad(slope)
  
  # Compute area of a cell (horizontal area or ground area)
  cell_horizontal_surface <- cell_size * cell_size
  cell_surface <- cell_horizontal_surface / cos(slope_rad)
  
  # Azimuth of the vector orthogonal to the ground in the x,y system
  bottom_azimut <- -aspect + north_to_x_cw
  bottom_azimut_rad <- deg2rad(bottom_azimut)
  
  # Prepare rays table
  rays_dt <- as.data.table(rays$rays)
  rays_dt[, `:=`(cos_heightangle = cos(height_angle),
                 sin_heightangle = sin(height_angle),
                 cos_azimut = cos(azimut),
                 sin_azimut = sin(azimut))]
  
  
  # For each ray, get relative cells that contains potential trees
  # (cells relative to any target cell enlighted by rays at its center)
  potcells_rays_rel <-
    sl_get_potentialcells_rays_relative(rays$rays,
                                        cell_size,
                                        slope_rad, bottom_azimut_rad,
                                        max_height = max(trees$height_m),
                                        max_cradius = max(trees$cradius_m))
  
  
  # Create table with all rays that enlight every potential cell center
  potcells_rays <- data.table(
    id_target = rep(cells$id_cell, each = nrow(potcells_rays_rel)),
    x_target = rep(cells$x_center, each = nrow(potcells_rays_rel)),
    y_target = rep(cells$y_center, each = nrow(potcells_rays_rel)),
    z_target = rep(cells$z_center, each = nrow(potcells_rays_rel)),
    
    id_ray = rep(potcells_rays_rel$id_ray, times = n_cells * n_cells),
    xid_rel = rep(potcells_rays_rel$xid_rel, times = n_cells * n_cells),
    yid_rel = rep(potcells_rays_rel$yid_rel, times = n_cells * n_cells)
  )
  
  
  # For each ray coming to a given target cell, get cells with associated shift
  # that contains trees that could potentially intercept the given ray
  potcells_rays[, c("xid_cell", "yid_cell",
                    "x_shift", "y_shift",
                    "outside") := cells_rel2abs(x_target, y_target,
                                                xid_rel, yid_rel,
                                                cell_size, n_cells)]
  
  potcells_rays[, z_shift := get_z(x_shift, y_shift, slope_rad, bottom_azimut_rad)]
  
  
  # Remove potential cells outside the plot if the torus system is disabled (shifted cells)
  if (!use_torus) {
    potcells_rays <- potcells_rays[outside == FALSE,]
  }
  
  
  # Add ray informations
  potcells_rays <- cells[potcells_rays, on = c(x_id = "xid_cell", y_id = "yid_cell")]
  potcells_rays <- rays_dt[potcells_rays, on = .(id_ray)]
  
  interceptions <- sl_get_interceptions(potcells_rays, trees)
  
  
  # METHOD TO OPIMIZE MEMORY WITH BASE R
  
  # Need to clusterize to avoid having very large dataframes (for memory)
  # cluster_size <- 100
  # id_targets <- 1:(n_cells*n_cells)
  # clusters <- split(id_targets, ceiling(seq_along(id_targets)/cluster_size))
  #
  # interceptions <- lapply(clusters,
  #                         function(c) {
  #                           sl_get_interceptions(potcells_rays[id_target %in% c],
  #                                                trees, rays_dt, use_torus)
  #                         })
  # interceptions <- rbindlist(interceptions)
  
  
  # Compute potential and intercepted energy for each interception
  if (turbid_medium) {
    interceptions <- sl_compute_energy_turbid(interceptions,
                                              slope_rad, bottom_azimut_rad,
                                              cell_surface)
  } else {
    interceptions <- sl_compute_energy_porous(interceptions,
                                              slope_rad, bottom_azimut_rad,
                                              cell_surface)
  }
  
  
  # Sum of diffuse and direct energy coming into slope surface
  e_slope_tot <- sum(rays$e_slope)
  
  # Summarize targets into two trees and cells datasets
  out <- sl_summarize_interceptions(interceptions, trees$id_tree,
                                    cells, cell_surface,
                                    e_slope_tot)
  
  return(out)
}
