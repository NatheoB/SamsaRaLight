#' Summarize interceptions output into trees and cells interception datasets
#'
#' @param out_targets Data.frame with bound targets output from sl_run_target()
#' @param trees_id vector of ids of the trees composing the stand
#' @param cells data.frame with id_cell and corresponding cell center coordinates
#' @param cell_surface Area of a cell considering slope (in m2)
#' @param e_slope_tot Sum of diffuse and direct energy coming into a cell considering slope (MJ)
#'
#' @return a list of trees and cells interception datasets
#'
#' @import data.table
#' @export
sl_summarize_interceptions <- function(interceptions, trees_id,
                                       cells,
                                       cell_surface_slope,
                                       e_slope_tot) {

  # Combine targets to have one global output per tree
  # Get mean energy for each tree
  out_trees <- interceptions[, .(epot = sum(epot),
                                 e = sum(e)), by = id_tree][, lci := (1 - e/epot)]

  # Energy of small trees (set lci to 1)
  # Small trees that have not intercepted energy (thus, not in interceptions dataset)
  interceptions_trees_id <- unique(interceptions$id_tree)
  if (length(interceptions_trees_id) < length(trees_id)) {

    smalltrees_id <- setdiff(trees_id, interceptions_trees_id)
    out_trees <- rbindlist(list(
      out_trees,
      data.table(
        id_tree = smalltrees_id,
        epot = 0,
        e = 0,
        lci = 1
      )
    ))

  }


  # Compute energy for each cell
  out_cells <- interceptions[, .(eintercepted_tot = sum(e)), by = id_target]

  out_cells[, e := e_slope_tot - (eintercepted_tot / cell_surface_slope)]
  out_cells[, erel := e / e_slope_tot]

  out_cells <- out_cells[cells, on = c(id_target = "id_cell")][, .(id_target, e, erel)]
  setnames(out_cells, "id_target", "id_cell")

  # Return both datasets as a list
  list("trees" = out_trees, "cells" = out_cells)
}


#' Compute potential and intercepted energy from interceptions (rays X trees)
#' using turbid medium (cf Ligot et al. 10.1139/cjfr-2013-0494)
#'
#' @import data.table
#' @export
sl_compute_energy_turbid <- function(interceptions,
                                     slope_rad, bottom_azimut_rad,
                                     cell_surface) {

  # Compute projection of energy on plane parallel to slope in MJ/m2.
  # And convert it into MJ per cell
  interceptions[, scalar_slope := cos(slope_rad) * sin(height_angle) +
                  sin(slope_rad) * cos(height_angle) * cos(azimut - bottom_azimut_rad)]
  interceptions[, e_incident_slope := scalar_slope * e_incident * cell_surface]

  # Compute potential energy from each ray intercepted trees
  interceptions[, epot := sl_beerlambert(incident_energy = e_incident_slope,
                                         path_length = length,
                                         leaf_area_density = crown_lad,
                                         extinction_coef = 0.5,
                                         clumping_factor = 1)]

  # Compute product of LAD and path length (used below)
  # (same to compute absorption of a ray with path length l in a volume of density d than
  # absorption of a ray with path length l*d in a volume of density 1)
  interceptions[, length_lad1 := length * crown_lad]

  # Arrange interception from the farthest to the nearest intercepted crown
  setorder(interceptions, -distance)

  # Compute total length of crown of volume 1 intercepted before hitting the given tree for each ray
  interceptions[, length_lad1_before := cumsum(length_lad1) - length_lad1, by = .(id_target, id_ray)]

  # Compute final intercepted energy

  # Energy of the ray that has been absorbed by previous trees when arriving to the given tree
  interceptions[, e_intercepted_before := sl_beerlambert(incident_energy = e_incident_slope,
                                                         path_length = length_lad1_before,
                                                         leaf_area_density = 1,
                                                         extinction_coef = 0.5,
                                                         clumping_factor = 1)]

  # Compute incident energy of a ray arriving to the tree after interception by all surrounding trees
  interceptions[, e_incident_after := e_incident_slope - e_intercepted_before]

  # Compute energy intercepted by the tree from a single ray after attenuation by previous trees
  interceptions[, e := sl_beerlambert(incident_energy = e_incident_after,
                                      path_length = length_lad1,
                                      leaf_area_density = 1,
                                      extinction_coef = 0.5,
                                      clumping_factor = 1)]

  # Remove trees outside the plot
  # Used to compute light attenuation when using torus system
  # interceptions <- interceptions[outside == FALSE, ]

  # Keep necessary columns
  interceptions[, .(id_target, id_ray, id_tree, epot, e)]
}


#' Compute potential and intercepted energy from interceptions (rays X trees)
#' using porous envelop (cf Ligot et al. 10.1139/cjfr-2013-0494)
#'
#' @import data.table
#' @export
sl_compute_energy_porous <- function(interceptions,
                                     slope_rad, bottom_azimut_rad,
                                     cell_surface) {

  # Compute projection of energy on plane parallel to slope in MJ/m2.
  # And convert it into MJ per cell
  interceptions[, scalar_slope := cos(slope_rad) * sin(height_angle) +
                  sin(slope_rad) * cos(height_angle) * cos(azimut - bottom_azimut_rad)]
  interceptions[, e_incident_slope := scalar_slope * e_incident * cell_surface]

  # Compute potential energy from each ray intercepted trees
  interceptions[, epot := (1 - crown_p) * e_incident_slope]

  # Arrange interception from the farthest to the nearest intercepted crown
  setorder(interceptions, -distance)

  # Compute total attenuation coeff before hitting the given tree for each ray
  interceptions[, total_attenuation := cumprod(crown_p) / crown_p, by = .(id_target, id_ray)]

  # Compute final intercepted energy

  # Compute incident energy of a ray arriving to the tree after interception by all surrounding trees
  interceptions[, e_incident_after := total_attenuation * e_incident_slope]

  # Compute energy intercepted by the tree from a single ray after attenuation by previous trees
  interceptions[, e := (1 - crown_p) * e_incident_after]

  # Remove trees outside the plot
  # Used to compute light attenuation when using torus system
  # interceptions <- interceptions[outside == FALSE, ]

  # Keep necessary columns
  interceptions[, .(id_target, id_ray, id_tree, epot, e)]
}


#' Search for interceptions between all rays coming to all target cells and trees
#'
#' @export
sl_get_interceptions <- function(potcells_rays, trees) {

  # Create the table with all possible trees intercepted by each ray (adding ray info)
  rays_trees <- potcells_rays[trees, on = .(id_cell), nomatch = NULL, allow.cartesian=TRUE]

  # Get characteristics of interception for each rayXtree (depending on crown form)
  # CAREFUL : 2nd shift is applied to shift coordinates to set target cell as origin

  # Ellipsoidal crowns
  rays_trees[crown_type == 2,
             c("length", "distance") :=
               sl_intercept_crown_ellipsoid(
                 cos_heightangle, sin_heightangle,
                 cos_azimut, sin_azimut,
                 x, y, z + cbh_m + (height_m - cbh_m)/2,
                 x_shift - x_target, y_shift - y_target, z_shift - z_target,
                 cradius_m, cradius_m, (height_m - cbh_m)/2)]

  # Paraboloidal crown
  rays_trees[crown_type == 1,
             c("length", "distance") :=
               sl_intercept_crown_paraboloid(
                 cos_heightangle, sin_heightangle,
                 cos_azimut, sin_azimut,
                 x, y, z + cbh_m,
                 x_shift - x_target, y_shift - y_target, z_shift - z_target,
                 cradius_m, cradius_m, height_m - cbh_m)]

  # Remove raysXtrees that have not intercept each other
  # And select necessary columns
  na.omit(rays_trees, cols="distance")[, .(id_target,
                                           id_ray, azimut, height_angle, e_incident,
                                           id_tree, crown_lad, crown_p,
                                           length, distance, outside)]
}


#' Get relative cell ids that could contain trees potentially intercepting a given ray
#'
#' @description The trees which can intercept a beam are located in cells with
#'  their center located in a competition rectangle of length L (beam X direction)
#'  + R (opposite direction) and width R (directions perpendicular X to beam).
#'
#' @param rays data.frame with n rows and 5 columns:
#' \itemize{
#'  \item{id}{Unique id of the ray}
#'  \item{azimut}{Azimut of the ray in radians}
#'  \item{height_angle}{Angle between beam and soil (in radians)}
#'  \item{e_init}{Initital energy of ray before crossing the canopy (in MJ.m-2)}
#'  \item{e_current}{Current energy of ray (in MJ.m-2)}
#'  \item{direct}{true if the ray is direct false if it is diffuse}
#' }
#' @param cell_size Length of the side of a squared cell composing the stand (in meters)
#' @param slope_rad Slope of the stand (in radians)
#' @param bottom_azimut_rad Azimut of the vector orthogonal to the ground in the x,y system (in radians)
#' @param max_height Max possible height of a tree in Samsara model (in meters)
#' @param max_cradius Max possible crown radius of a tree in Samsara model (in meters)
#'
#' @return Return a data.frame with n rows for each ray X possible cell and 3 columns
#' \itemize{
#'  \item{id_ray}{Unique id of the ray}
#'  \item{xid_rel}{Relative X-axis id compared to the target cell (cell where ray is targeted to its center)}
#'  \item{yid_rel}{Relative Y-axis id compared to the target cell (cell where ray is targeted to its center)}
#' }
#'
#' @export
sl_get_potentialcells_rays_relative <- function(rays,
                                                cell_size,
                                                slope_rad,
                                                bottom_azimut_rad,
                                                max_height,
                                                max_cradius) {
  rays %>%
    dplyr::mutate(

      # Computes lateral = the boundary to add to the competition
      # rectangle to take into account cells center
      # instead of trees position.
      # The boundary depends on beam azimut.
      azt = case_when(
        azimut < pi / 4 ~ azimut,
        (azimut >= pi / 4) & (azimut < pi / 2) ~ pi / 2 - azimut,
        (azimut >= pi / 2) & (azimut < 3 * pi / 4) ~ azimut - pi / 2,
        (azimut >= 3 * pi / 4) & (azimut < pi) ~ pi - azimut,
        (azimut >= pi) & (azimut < 5 * pi / 4) ~ azimut - pi,
        (azimut >= 5 * pi / 4) & (azimut < 3 * pi / 2) ~ 3 * pi / 2 - azimut,
        (azimut >= 3 * pi / 2) & (azimut < 7 * pi / 4) ~ azimut - 3 * pi / 2,
        azimut >= 7 * pi / 4 ~ 2 * pi - azimut,
        TRUE ~ 0
      ),
      lateral = cell_size / sqrt(2) * sin(azt + pi / 4),

      # Beam width = max lateral distance from the beam to a cell center
      # able to own a tree which can intercept the beam.
      R = max_cradius + lateral,

      # Beam reach maximum distance along the beam beyond which the cells
      # cannot own trees which can intercept the beam (too high).
      L = max_height / (tan(height_angle) + cos(azimut - bottom_azimut_rad) * tan(slope_rad)) + lateral,

      # Coordinates of the four corners of the competition rectangle.
      sinA = sin(azimut),
      cosA = cos(azimut),

      x1 = R * sinA + L * cosA,
      y1 = L * sinA - R * cosA,
      x2 = L * cosA - R * sinA,
      y2 = L * sinA + R * cosA,
      x3 = R * (sinA - cosA),
      y3 = -R * (sinA + cosA),
      x4 = -R * (sinA + cosA),
      y4 = R * (cosA - sinA),

      x_min = pmin(x1, x2, x3, x4),
      x_max = pmax(x1, x2, x3, x4),
      y_min = pmin(y1, y2, y3, y4),
      y_max = pmax(y1, y2, y3, y4),

      # Round into relative-to-target coordinates of the cell center in which x/y_min/max are
      x_min = ceiling(x_min / cell_size) * cell_size,
      x_max = floor(x_max / cell_size) * cell_size,
      y_min = ceiling(y_min / cell_size) * cell_size,
      y_max = floor(y_max / cell_size) * cell_size,

      # Number of cells between min and max for both axis x and y
      nx = floor((x_max - x_min) / cell_size + 0.5) + 1,
      ny = floor((y_max - y_min) / cell_size + 0.5) + 1
    ) %>%

    # Add a row for each candidate cell for each ray
    tidyr::uncount(nx, .remove = F) %>%
    tidyr::uncount(ny, .remove = F) %>%
    dplyr::group_by(id_ray) %>%
    dplyr::mutate(
      # Id of the cell from the minimum ones
      xid = rep(1:unique(nx), times = unique(ny)) - 1,
      yid = rep(1:unique(ny), each = unique(nx)) - 1,

      # Coordinates of the candidate cells (be careful: y coordinates is inverse direction of y grid)
      # id = 1 is the greater y coordinates
      x = x_min + xid * cell_size,
      y = y_max - yid * cell_size,

      # Relative id of the cell compared to the origin
      xid_rel = floor(x / cell_size),
      yid_rel = - floor(y / cell_size)

    ) %>%
    dplyr::ungroup() %>%

    # Check if candidate cell is located inside the competition rectangle
    dplyr::mutate(
      is_inside = (x * sinA - y * cosA < R) &
        (x * cosA + y * sinA < L) &
        (-x * sinA + y * cosA < R) &
        (x * cosA + y * sinA > -R)
    ) %>%

    # Return only raysXcells that are possible
    dplyr::filter(is_inside) %>%
    dplyr::select(id_ray, xid_rel, yid_rel)
}



#' Compute intercepted energy using Beer Lambert law
#'
#' @param incident_energy energy that arrive to the crown
#' @param path_length length of the path within the crown (in meters)
#' @param leaf_area_density LAD of the crown (m2 of leaves per m3 of crown volume)
#' @param extinction_coef Probability of a leaf to intercept the ray (linked to leaf orientation)
#' @param clumping_factor Aggregation of leaves within the crown volume (1 is homogeneous)
#'
sl_beerlambert <- function(incident_energy,
                           path_length,
                           leaf_area_density = 0.5,
                           extinction_coef = 0.5,
                           clumping_factor = 1) {

  incident_energy * (1 - exp(- extinction_coef * clumping_factor * leaf_area_density * path_length))

}
