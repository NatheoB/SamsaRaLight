#' Center the surveyed trees of the core polygon in a rectangle plot
#'  
#' @param trees a data.frame with one row for each tree and their position in a (x,y) grid.
#'  It must be composed of two columns "x" and "y"
#' @param cell_size double - Length of the side of a squared cell composing the stand (in meters)
#' @param core_polygon_df data.frame - Coordinates of the vertices composing the core polygon.
#'  It must be composed of two columns "x" and "y" 
#' @param use_rect_zone boolean - How to consider the delimitation of the inventory zone,
#'  used to define the square plot and when fill_around is TRUE, for filling the trees aorund.
#'  Either to use the given or computed core polygon (FALSE), 
#'  or the minimum-area enclosing rectangle from the polygon (TRUE).
#' @param fill_around boolean - If TRUE, add trees around the core polygon within the stand 
#'  until the stand reach the same total basal area per hectare than the core polygon
#'  
#'
#' @import sf sfheaders concaveman lwgeom
#'
#' @export
#'
create_rect_stand <- function(trees, 
                              cell_size,
                              core_polygon_df = NULL,
                              use_rect_zone = FALSE,
                              fill_around = FALSE) {
  
  
  # DEFINE THE CORE INVENTORY ZONES ----
  
  ## Check if the core polygone is a polygone (at least 3 tops) ----
  if (!is.null(core_polygon_df) && nrow(core_polygon_df) < 3) {
    message("WARNING: the core polygone has less than 3 tops: we did not considered it")
    core_polygon_df <- NULL
  }
  
  
  ## Define a polygon that encompasses all the trees ----
  # Either the polygon is initially given, or created from concaveman if it was NULL
  coords_sf <- st_as_sf(trees, coords = c("x", "y")) 
  
  if (is.null(core_polygon_df)) {
    core_polygon_sf <- concaveman::concaveman(coords_sf)[[1]]
    core_polygon_df <- st_coordinates(core_polygon_sf) %>% 
      as.data.frame() %>% 
      dplyr::select(x = X, y = Y)
  } else {
    core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  }
  
  
  ## Ensure that all trees are in the core polygon ----
    # use st_intersect and not st_within to also consider trees in the edges of the polygon
  if (!all(st_intersects(coords_sf, core_polygon_sf, sparse = FALSE))) {
    stop("ERROR: some trees are outside the core polygon")
  }
  
  
  ## Get the minimum-area enclosing rectangle from the polygon ----
  rotated_rect_sf <- st_minimum_rotated_rectangle(core_polygon_sf)
  
  rotated_rect_df <- st_coordinates(rotated_rect_sf) %>% 
    as.data.frame() %>% 
    dplyr::select(x = X, y = Y) %>% 
    dplyr::distinct() # Because 5 vertices (the last one is the same as the first one)
  
    # To update the sf object vertices after having removed the last vertex
  rotated_rect_sf <- sfheaders::sf_polygon(rotated_rect_df)
  
  
  ## Set the main inventory zone ----
  if (use_rect_zone) {
    inv_zone_df <- rotated_rect_df
    inv_zone_sf <- rotated_rect_sf
  } else {
    inv_zone_df <- core_polygon_df
    inv_zone_sf <- core_polygon_sf
  }
  
  
  

  # CREATE THE SQUARE PLOT ----
  
  ## Get the range of the inventory zone ----
  inv_size_x <- max(inv_zone_df$x) - min(inv_zone_df$x)
  inv_size_y <- max(inv_zone_df$y) - min(inv_zone_df$y)
  
  
  ## Find the minimum of cells in the X- and Y-axis ----
  n_cells_x <- inv_size_x %/% cell_size + 1
  n_cells_y <- inv_size_y %/% cell_size + 1
  
  
  ## Compute the square plot size ----
  plot_size_x <- cell_size * n_cells_x
  plot_size_y <- cell_size * n_cells_y
  
  
  ## Find the shift to apply ---
  # on base coordinates to center the inventory within the square plot
  diff_x <- plot_size_x - inv_size_x
  diff_y <- plot_size_y - inv_size_y
  
  shift_x <- - ( min(inv_zone_df$x) - diff_x / 2 )
  shift_y <- - ( min(inv_zone_df$y) - diff_y / 2 )
  
  
  ## Shift coordinates ----
  trees_shifted <- trees %>% 
    dplyr::mutate(
      x = x + shift_x,
      y = y + shift_y
    ) %>% 
    as.data.frame()
  
  inv_zone_shifted_df <- inv_zone_df %>% 
    dplyr::mutate(
      x = x + shift_x,
      y = y + shift_y
    ) %>% 
    as.data.frame()
  
  inv_zone_shifted_sf <- sfheaders::sf_polygon(inv_zone_shifted_df)
  
    
  
  
  # FILL AROUND THE CORE INVENTORY ZONE ----
  
  ## Compute total basal area per hectare in the polygon ----
  core_area_m2 <- st_area(inv_zone_shifted_sf)
  core_area_ha <- core_area_m2 / 10000
  
  batot_target_m2ha <- sum( 3.1416 * (trees_shifted$dbh_cm / 200) ^ 2 ) / core_area_ha
  
  
  ## Fill with trees outside the polygon ----
  plot_size_x <- cell_size * n_cells_x
  plot_size_y <- cell_size * n_cells_y
  
  plot_area_ha <- (plot_size_x * plot_size_y) / 10000
  
  current_batot_m2ha <- sum( 3.1416 * (trees_shifted$dbh_cm / 200) ^ 2 ) / plot_area_ha
  new_id <- max(trees_shifted$id_tree) + 1
  
  i <- 0
  i_max <- 1e6 # Maximum iterations 
  
  trees_to_add <- list() # Dataframe of new trees to add withhin the rectangle
  
  if (fill_around) {
    
    # until we reach the target basal area
    while ((current_batot_m2ha < batot_target_m2ha) & (i < i_max)) {
      
      # print(i)
      
      # Get a random tree
      tree_to_add <- trees_shifted[sample(nrow(trees_shifted), 1), ]
      
      # Create a new id
      tree_to_add$id_tree <- new_id
      new_id <- new_id + 1
      
      # Add a random position until it is outside
      new_x <- runif(1, min = 0, max = plot_size_x)
      new_y <- runif(1, min = 0, max = plot_size_y)
      while (st_intersects(st_point(c(new_x, new_y)), 
                           inv_zone_shifted_sf, sparse = FALSE)[1,1]) {
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
    
  }
  
  trees_to_add <- dplyr::bind_rows(trees_to_add)
  
  
  ## Final inventory tree ----
  # Specify that the trees added are not an original ones (not inside the core polygon)
  trees_filled <- dplyr::bind_rows(
    trees_shifted %>% dplyr::mutate(added_to_fill = FALSE), 
    trees_to_add %>% dplyr::mutate(added_to_fill = TRUE)
  ) %>% 
    as.data.frame()
  
  
  # OUTPUT A LIST  ----
  # about the new dataset, the core polygon and the plot geometry
  
  return(list(
    "trees" = trees_filled,
    "inv_zone_df" = inv_zone_shifted_df,
    "inv_zone_sf" = inv_zone_shifted_sf,
    "info" = list("core_area_ha" = core_area_ha, 
                  "core_batot_m2ha" = batot_target_m2ha,
                  "n_added_tree" = i,
                  "new_area_ha" = plot_area_ha,
                  "new_batot_m2ha" = current_batot_m2ha,
                  "cell_size" = cell_size,
                  "n_cells_x" = n_cells_x, 
                  "n_cells_y" = n_cells_y,
                  "shift_x" = shift_x,
                  "shift_y" = shift_y)
  ))
}