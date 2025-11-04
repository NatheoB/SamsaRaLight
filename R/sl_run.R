#' Compute samsara ligth radiative balance
#'
#' @param trees a data.frame with one row for each tree and 14 columns describing tree composing the stand:
#' \itemize{
#'  \item{id: }{Unique id of the tree}
#'  \item{x: }{X position of the tree within the stand}
#'  \item{y: }{Y position of the tree within the stand}
#'  \item{dbh_cm: }{Diameter at breast height (1.3m) of the tree trunk (in cm)}
#'  \item{crown_type: }{
#'    \itemize{
#'      \item{"E": symetric ellispoidal crown with the maximum crown radius at the middle of the crown depth.
#'        Crown radius is the mean value of the radius towards the fourth cardinal points}
#'      \item{"P": symetric paraboloidal crown with the maximum crown radius at the crown base.
#'        crown radius is the mean value of the radius towards the fourth cardinal points}  
#'      \item{"2E": ellipsoidal crown composed of two semi-ellipsoid (one above and one below).
#'        Same as "E" but height of each semi-ellipsoid is given by hmax_m variable}
#'      \item{"8E": ellispoidal crown composed of eight eighth of ellipsoids (four above and four below). 
#'        Same as "2E" but radius of each eight of ellipsoid is given by rn_m, rs_m, re_m, rw_m variables.}
#'      \item{"4P": paraboloidal crown composed of four fourth of paraboloid.
#'        Maximum crown radius at the crown base and radius of each fourth of paraboloid is given by rn_m, rs_m, re_m, rw_m variables.}  
#'    }
#'  }
#'  \item{h_m: }{Height of the tree (in meters)}
#'  \item{hbase_m: }{Height of the base of the tree crown (in meters)}
#'  \item{hmax_m: }{Height of the maximum radius of the tree crown (in meters)}
#'  \item{rn_m: }{Radius of the tree crown toward Northern direction (in meters)}
#'  \item{rs_m: }{Radius of the tree crown toward Southern direction (in meters)}
#'  \item{re_m: }{Radius of the tree crown toward Eastern direction (in meters)}
#'  \item{rw_m: }{Radius of the tree crown toward Western direction (in meters)}
#'  \item{crown_openess: }{Crown openess of the crown (no unit) when considering porous envelop 
#'    (i.e. constant proportion of light energy interception when a ray intercepts a crown)}
#'  \item{crown_lad: }{Leaf Area Density (in m2 of leaves per m3 of crown) when considering turbid medium 
#'    (i.e. density parameter of Beer Lambert law when a ray intercepts a crown)}
#' }
#' @param monthly_rad data.frame - Monthly horizontal radiation (Hrad) and diffuse to global ratio (DGratio)
#'    Computed with function `samsaRaLight::sl_get_monthlyrad()`
#' @param sensors data.frame - Position of the sensor within the stand. Can be set to NULL if no sensor.
#' \itemize{
#'  \item{id: }{Unique id of the sensor}
#'  \item{x: }{X position of the sensor within the stand}
#'  \item{y: }{Y position of the sensor within the stand}
#'  \item{h_m: }{Height above ground of the sensor (in meters)}
#' }
#' @param sensors_only boolean - To compute interception only for the sensors ?
#' @param latitude double - Latitude of the plot (in degrees)
#' @param start_day integer between 1 and 365 - First day of the vegetative period
#' @param end_day integer between 1 and 365 - Last day of the vegetative period
#' @param slope double - Slope of the plot (in degrees)
#' @param north_to_x_cw double - Angle from North to x axis clockwise. (in degrees)
#'    Default correspond to a Y axis oriented toward the North.
#' @param aspect double - Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
#'    northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
#' @param cell_size double - Length of the side of a squared cell composing the stand (in meters)
#' @param n_cells_x integer - Number of cells of the X-axis side of the rectangle stand
#' @param n_cells_y integer - Number of cells of the Y-axis side of the rectangle stand
#' @param soc boolean - Standard Overcast Sky, if false: Uniform Overcast Sky
#' @param height_anglemin double - Angle minimum between beam and soil (in degrees)
#' @param direct_startoffset double - Angle at which to start first direct ray (in degrees)
#' @param direct_anglestep double - Hour angle between two direct beams (in degrees)
#' @param diffuse_anglestep double - Hour angle between two diffuse beams (in degrees)
#' @param use_torus if True, compute light competition using borders modelled with a
#'  torus system, otherwise,borders are open grasslands
#' @param turbid_medium If TRUE, crown are considered as turbid medium, otherwise,
#'  considered as porous envelope
#' @param trunk_interception Consider interception of rays by trunks
#' @param detailed_output boolean - If TRUE, sensors/cells/trees outputs also contain diffuse/direct energies, 
#'  and sensors/cells outputs contain energies both on the slope and on a horizontal plane.
#'  If FALSE, energy given is only the total energy (sum of diffuse and direct) on the slope for trees and cells
#'  and on a horizontal plane for sensors.
#'
#' @import data.table
#' @importFrom Rcpp sourceCpp
#' @useDynLib SamsaRaLight, .registration = TRUE
#'
#' @export
#'
sl_run <- function(trees,
                   monthly_rad,
                   sensors = NULL, sensors_only = FALSE,
                   latitude = 46,
                   start_day = 1, end_day = 365,
                   slope = 0,
                   north_to_x_cw = 90,
                   aspect = 0,
                   cell_size = 10,
                   n_cells_x = 10,
                   n_cells_y = 10,
                   soc = TRUE,
                   height_anglemin = 10,
                   direct_startoffset = 0,
                   direct_anglestep = 5,
                   diffuse_anglestep = 15,
                   use_torus = TRUE,
                   turbid_medium = TRUE,
                   trunk_interception = TRUE,
                   detailed_output = FALSE) {

  # Checks the input arguments
  check_coordinates(trees, sensors, cell_size, n_cells_x, n_cells_y)
  
  
  # Create monthly rays
  monthly_rays <- sl_create_monthly_rays(monthly_rad = monthly_rad,
                                         latitude = latitude,
                                         start_day = start_day,
                                         end_day = end_day,
                                         soc = soc,
                                         slope = slope,
                                         north_to_x_cw = north_to_x_cw,
                                         aspect = aspect,
                                         height_anglemin = height_anglemin,
                                         direct_startoffset = direct_startoffset,
                                         direct_anglestep = direct_anglestep,
                                         diffuse_anglestep = diffuse_anglestep)
  
  # Run call to c++ script
  out <- sl_run_rcpp(
    trees, 
    sensors, sensors_only,
    monthly_rays$rays, 
    monthly_rays$energies[["slope_direct"]], 
    monthly_rays$energies[["slope_diffuse"]],
    monthly_rays$energies[["horizontal_direct"]], 
    monthly_rays$energies[["horizontal_diffuse"]],
    slope, north_to_x_cw, aspect,
    cell_size, n_cells_x, n_cells_y,
    use_torus, turbid_medium, trunk_interception)

  
  # Filter the output datasets if asked
  if (!detailed_output) {
    
    # Sensors dataset
    out$sensors <- out$sensors %>% 
      dplyr::select(id_sensor, x, y, z, 
                    e = e_horizontal, pacl = pacl_horizontal, punobs = punobs_horizontal) %>% 
      as.data.frame()
    
    # Cells dataset
    out$cells <- out$cells %>% 
      dplyr::select(id_cell, x_center, y_center, z_center, 
                    e = e_slope, pacl = pacl_slope, punobs = punobs_slope) %>% 
      as.data.frame()
    
    # Trees dataset
    out$trees <- out$trees %>% 
      dplyr::select(id_tree, x, y, z,
                    epot, e, lci, eunobs) %>% 
      as.data.frame()
    
    # Remove the interception matrices
    out$interceptions <- NULL
  }
  
  # Create the SamsaRaLight output object
  list(
    "input" = list(
      "trees" = trees,
      "monthly_rad" = monthly_rad,
      "info" = list(
        "latitude" = latitude,
        "start_day" = start_day, 
        "end_day" = end_day,
        "slope" = slope,
        "north_to_x_cw" = north_to_x_cw,
        "aspect" = aspect,
        "cell_size" = cell_size,
        "n_cells_x" = n_cells_x,
        "n_cells_y" = n_cells_y,
        "soc" = soc,
        "height_anglemin" = height_anglemin,
        "direct_startoffset" = direct_startoffset,
        "direct_anglestep" = direct_anglestep,
        "diffuse_anglestep" = diffuse_anglestep,
        "use_torus" = use_torus,
        "turbid_medium" = turbid_medium,
        "trunk_interception" = trunk_interception
      )
    ),
    "monthly_rays" = monthly_rays,
    "output" = out
  )
}



#' Check if all the trees and sensors are within the plot limits
#' @noRd
check_coordinates <- function(trees, sensors, cell_size, n_cells_x, n_cells_y) {
  
  # Check trees coordinates
  outside_trees <- trees %>% 
    dplyr::mutate(
      is_outside = x < 0 | y < 0 | x > cell_size * n_cells_x | y > cell_size * n_cells_y
    ) %>% 
    dplyr::pull(is_outside)
  
  if (sum(outside_trees) > 0) stop("Some trees are outside the plot limits...")
  
  
  # Check sensors coordinates
  if (!is.null(sensors)) {
    outside_sensors <- sensors %>% 
      dplyr::mutate(
        is_outside = x < 0 | y < 0 | x > cell_size * n_cells_x | y > cell_size * n_cells_y
      ) %>% 
      dplyr::pull(is_outside)
    
    if (sum(outside_sensors) > 0) stop("Some sensors are outside the plot limits...")
  }
  
}
