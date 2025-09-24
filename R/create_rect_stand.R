#' Center the surveyed trees of the core polygon in a rectangle plot
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
#' @param cell_size double - Length of the side of a squared cell composing the stand (in meters)
#' @param core_polygon data.frame - Coordinates of the ertices composing the core polygon.
#'  It must be composed of two columns "x" and "y" 
#'
#' @import sf sfheaders concaveman
#'
#' @export
#'
create_rect_stand <- function(trees, cell_size, core_polygon = NULL) {
  
  # Check if the core polygone is a polygone (at least 3 tops)
  if (!is.null(core_polygon) && nrow(core_polygon) < 3) {
    message("The core polygone has less than 3 tops: we did not considered it")
    core_polygon <- NULL
  }
  
  # Get the square size ----
  if (is.null(core_polygon)) {
    
    # Take the range of surveyed trees
    inv_size_x <- max(trees$x) - min(trees$x)
    inv_size_y <- max(trees$y) - min(trees$y)
    
  } else {
    
    # Take range of core polygon if it is not null
    inv_size_x <- max(core_polygon$x) - min(core_polygon$x)
    inv_size_y <- max(core_polygon$y) - min(core_polygon$y)
    
  }
  
  
  # Find the minimum of cells in the X- and Y-axis ----
  n_cells_x <- inv_size_x %/% cell_size + 1
  n_cells_y <- inv_size_y %/% cell_size + 1
  
  
  # Compute the square plot size ----
  plot_size_x <- cell_size * n_cells_x
  plot_size_y <- cell_size * n_cells_y
  
  
  # Find the shift to apply ---
  # on base coordinates to center the inventory within the square plot
  diff_x <- plot_size_x - inv_size_x
  diff_y <- plot_size_y - inv_size_y
  
  if (is.null(core_polygon)) {
    
    # Shift is based on the surveyed trees limits
    shift_x <- - ( min(trees$x) - diff_x / 2 )
    shift_y <- - ( min(trees$y) - diff_y / 2 )
  
  } else {
    
    # Otherwise, take the core polygone limits if it is defined
    shift_x <- - ( min(core_polygon$x) - diff_x / 2 )
    shift_y <- - ( min(core_polygon$y) - diff_y / 2 )
    
  }
  
  
  # Shift coordinates ----
  trees_shifted <- trees %>% 
    dplyr::mutate(
      x = x + shift_x,
      y = y + shift_y
    ) %>% 
    as.data.frame()
  
  core_polygon_shifted <- core_polygon
  if (!is.null(core_polygon_shifted)) {
    core_polygon_shifted <- core_polygon_shifted %>% 
      dplyr::mutate(
        x = x + shift_x,
        y = y + shift_y
      ) %>% 
      as.data.frame()
  }
  
  
  # Ensure that all trees are in the squared plot ----
  # It is not ensuired if the trees are iniially outside the given core polygon
  if (min(trees_shifted$x) < 0 || max(trees_shifted$x) > plot_size_x ||
      min(trees_shifted$y) < 0 || max(trees_shifted$y) > plot_size_y) {
    
    stop("Some trees are outside the core polygon")
    
  }
  
  # Return the plot, with some information
  list(
    "new_inventory" = trees_shifted,
    "core_polygon" = core_polygon_shifted,
    "info" = list("cell_size" = cell_size,
                  "n_cells_x" = n_cells_x, 
                  "n_cells_y" = n_cells_y,
                  "shift_x" = shift_x,
                  "shift_y" = shift_y)
  )
}