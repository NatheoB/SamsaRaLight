#' Create an equivalent rectangle plot around a core polygon
#'  composed only of the tree inventory.
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
#' @param n_cells_x integer - Number of cells of the X-axis side of the rectangle stand
#' @param n_cells_y integer - Number of cells of the Y-axis side of the rectangle stand
#' @param core_polygon data.frame - Coordinates of the ertices composing the core polygon.
#'  It must be composed of two columns "x" and "y" 
#'
#' @import sf sfheaders concaveman
#'
#' @export
#'
fill_square_inventory <- function(trees, 
                                  cell_size,
                                  n_cells_x, 
                                  n_cells_y,
                                  core_polygon = NULL)
  
{
  # Find a polygon that encompasses all the trees
  if (is.null(core_polygon)) {
    coords_sf <- st_as_sf(trees, coords = c("x", "y")) 
    core_polygon <- concaveman::concaveman(coords_sf)[[1]]
  } else {
    core_polygon <- sfheaders::sf_polygon(core_polygon)
  }
  
  # Compute total basal area per hectare in the polygon
  core_area_m2 <- st_area(core_polygon)
  core_area_ha <- core_area_m2 / 10000
  
  batot_target_m2ha <- sum( 3.1416 * (trees$dbh_cm / 200) ^ 2 ) / core_area_ha

  
  # Add trees outside the polygon, until we reach the target basal area
  plot_size_x <- cell_size * n_cells_x
  plot_size_y <- cell_size * n_cells_y
  
  plot_area_ha <- (plot_size_x * plot_size_y) / 10000
  
  current_batot_m2ha <- sum( 3.1416 * (data_trees$dbh_cm / 200) ^ 2 ) / plot_area_ha
  new_id <- max(data_trees$id_tree) + 1
  
  i <- 0
  i_max <- 1e6 # Maximum iterations 
  
  trees_to_add <- list() # Dataframe of new trees to add withhin the rectangle
  while ((current_batot_m2ha < batot_target_m2ha) & (i < i_max)) {
    
    # print(i)
    
    # Get a random tree
    tree_to_add <- data_trees[sample(nrow(data_trees), 1), ]
    
    # Create a new id
    tree_to_add$id_tree <- new_id
    new_id <- new_id + 1
    
    # Add a random position until it is outside
    new_x <- runif(1, min = 0, max = plot_size_x)
    new_y <- runif(1, min = 0, max = plot_size_y)
    while (st_within(st_point(c(new_x, new_y)), core_polygon, sparse = FALSE)[1,1]) {
      new_x <- runif(1, min = 0, max = plot_size_x)
      new_y <- runif(1, min = 0, max = plot_size_y)
    }
    tree_to_add$x <- new_x
    tree_to_add$y <- new_y
    
    # Add the tree to the dataset
    trees_to_add[[i+1]] <- tree_to_add
    
    # Increment the total basal area and the loop index
    current_batot_m2ha <- current_batot_m2ha + 3.1416 * (tree_to_add$dbh_cm / 200) ^ 2 / plot_area_ha
    i <- i + 1
  }
  trees_to_add <- dplyr::bind_rows(trees_to_add)
  

  # Final inventory tree
  # Specify that the trees added are not an original ones (not inside the core polygon)
  trees_filled <- dplyr::bind_rows(
    trees %>% dplyr::mutate(added_to_fill = FALSE), 
    trees_to_add %>% dplyr::mutate(added_to_fill = TRUE)
  ) %>% 
    as.data.frame()
  
  # Output a list with information about the new dataset and the core polygon
  return(list(
    "new_inventory" = trees_filled,
    "core_polygon" = core_polygon,
    "info" = list("core_area_ha" = core_area_ha, 
                  "core_batot_m2ha" = batot_target_m2ha,
                  "n_added_tree" = i,
                  "new_batot_m2ha" = current_batot_m2ha,
                  "cell_size" = cell_size,
                  "n_cells_x" = n_cells_x, 
                  "n_cells_y" = n_cells_y)
  ))
}