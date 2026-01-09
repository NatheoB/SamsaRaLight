#' Create a rectangular virtual stand from a tree inventory
#'
#' This function builds a rectangular (or square) virtual forest stand from a
#' user-provided tree inventory. Trees are spatially shifted so that the
#' inventory zone is centered within a regular grid of square cells.
#' Optionally, additional trees can be added around the core inventory area
#' to match its basal area per hectare.
#'
#' The function supports sloping terrain and coordinate system rotation, and
#' returns a fully prepared stand ready for use in the ray-tracing pipeline
#' (see \link{run_sl}).
#'
#' @param trees_inv A data.frame with one row per tree.
#'   See \link{check_inventory} for the required structure and validated columns.
#' @param cell_size Numeric. Side length of square cells composing the stand (meters).
#' @param slope Numeric. Slope of the plot (degrees).
#' @param aspect Numeric. Aspect of the slope, defined as the azimuth of the
#'   downslope direction, clockwise from North (degrees).
#'   North = 0, East = 90, South = 180, West = 270.
#' @param north2x Numeric. Clockwise angle from North to the X-axis (degrees).
#'   A value of 0 corresponds to a Y-axis oriented toward North.
#' @param sensors Optional data.frame defining position and height of the sensor within the stand.
#'   See \link{check_sensors} for the required structure and validated columns.
#' @param core_polygon_df Optional data.frame defining the core inventory polygon.
#'   Must contain columns \code{x} and \code{y}. If \code{NULL}, a concave hull
#'   is automatically computed from tree positions.
#' @param use_rect_zone Logical. If \code{TRUE}, the inventory zone is defined
#'   by the minimum-area enclosing rectangle of the core polygon.
#'   If \code{FALSE}, the core polygon itself is used.
#' @param fill_around Logical. If \code{TRUE}, trees are added outside the core
#'   polygon until the basal area per hectare of the full stand matches that
#'   of the core inventory.
#' @param verbose Logical. If \code{TRUE} (default), messages and warnings are
#'   printed during processing. If \code{FALSE}, output is silent.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{trees}}{
#'     Data.frame of the final tree inventory, including added trees if
#'     \code{fill_around = TRUE}.  
#'     The structure matches the inventory format validated by
#'     \link{check_inventory}, with additional derived variables required
#'     for ray tracing. A logical column \code{added_to_fill} indicates whether
#'     each tree originates from the initial inventory or was added to fill
#'     around the inventory zone.
#'   }
#'   \item{\code{core_polygon}}{
#'     List describing the inventory zone:
#'     \itemize{
#'       \item \code{df}: data.frame of polygon vertices
#'       \item \code{sf}: corresponding \code{sf} POLYGON
#'       \item \code{use_rect_zone}: logical flag
#'     }
#'   }
#'   \item{\code{transform}}{
#'     List of transformation and filling information, including core area,
#'     target and final basal area, number of added trees, and applied spatial shifts.
#'   }
#'   \item{\code{geometry}}{
#'     List describing stand geometry and terrain parameters
#'     (cell size, number of cells, slope, aspect, and orientation).
#'   }
#' }
#'
#' @details
#' The returned \code{trees} data.frame conforms to the inventory format checked
#' by \link{check_inventory}, with the following controlled modifications:
#' \itemize{
#'   \item Tree vertical position \code{z} is computed from terrain slope and aspect.
#'   \item Crown maximum radius height \code{hmax_m} is computed when fixed by crown
#'     geometry conventions:
#'     \itemize{
#'       \item \code{"P"} and \code{"4P"}: \code{hmax_m = hbase_m}
#'       \item \code{"E"} and \code{"4E"}: \code{hmax_m = hbase_m + 0.5 * (h_m - hbase_m)}
#'     }
#'     For crown types \code{"2E"} and \code{"8E"}, \code{hmax_m} must be provided by the user.
#'   \item If missing, column \code{dbh_cm} is added and filled with \code{NA}.
#'   \item If missing, crown interception properties (e.g. \code{crown_openness},
#'     \code{crown_lad}) are added using default values.
#' }
#'
#' The function ensures that all trees fall within the core inventory polygon,
#' applying small buffers if necessary to handle numerical precision issues.
#' Invalid polygons are automatically repaired when possible.
#'
#' When \code{fill_around = TRUE}, trees are randomly sampled from the original
#' inventory and positioned outside the core polygon until the target basal area
#' per hectare is reached for the full rectangular stand.
#'
#' @importFrom sf st_as_sf st_is_valid st_make_valid st_area st_buffer st_intersects st_minimum_rotated_rectangle
#' @importFrom sf st_point st_coordinates st_collection_extract
#' @importFrom concaveman concaveman
#' @importFrom sfheaders sf_polygon
#' @importFrom dplyr bind_rows mutate select distinct arrange
#' @importFrom tidyr expand_grid
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' data_prenovel <- SamsaRaLight::data_prenovel
#' trees_inv <- data_prenovel$trees
#'
#' stand <- create_sl_stand(
#'   trees_inv = trees_inv,
#'   cell_size = 5,
#'   slope = 10,
#'   aspect = 180,
#'   north2x = 0,
#'   use_rect_zone = TRUE,
#'   fill_around = FALSE,
#'   verbose = TRUE
#' )
#'
#' head(stand$trees)
#' }
#'
#' @export
create_sl_stand <- function(trees_inv, 
                            cell_size,
                            slope, 
                            aspect, 
                            north2x,
                            sensors = NULL,
                            core_polygon_df = NULL,
                            use_rect_zone = FALSE,
                            fill_around = FALSE,
                            verbose = TRUE) {
  
  # ARGUMENT CHECKS ----
  
  # trees_inv
  if (!inherits(trees_inv, "data.frame")) {
    stop("`trees_inv` must be a data.frame verified by check_inventory().", call. = FALSE)
  }
  
  # cell_size
  if (!is.numeric(cell_size) || length(cell_size) != 1 || is.na(cell_size)) {
    stop("`cell_size` must be a single numeric value.", call. = FALSE)
  }
  if (cell_size <= 0 || cell_size != as.integer(cell_size)) {
    stop("`cell_size` must be a positive integer (meters).", call. = FALSE)
  }
  
  # slope
  if (!is.numeric(slope) || length(slope) != 1 || is.na(slope)) {
    stop("`slope` must be a single numeric value (degrees).", call. = FALSE)
  }
  if (slope < 0 || slope >= 90) {
    stop("`slope` must be between 0 (inclusive) and 90 (exclusive) degrees.", call. = FALSE)
  }
  
  # aspect
  if (!is.numeric(aspect) || length(aspect) != 1 || is.na(aspect)) {
    stop("`aspect` must be a single numeric value (degrees).", call. = FALSE)
  }
  if (aspect < 0 || aspect >= 360) {
    stop("`aspect` must be between 0 (inclusive) and 360 (exclusive) degrees.", call. = FALSE)
  }
  
  # north2x
  if (!is.numeric(north2x) || length(north2x) != 1 || is.na(north2x)) {
    stop("`north2x` must be a single numeric value (degrees).", call. = FALSE)
  }
  if (north2x < 0 || north2x >= 360) {
    stop("`north2x` must be between 0 (inclusive) and 360 (exclusive) degrees.", call. = FALSE)
  }
  
  # core_polygon_df
  if (!is.null(core_polygon_df)) {
    if (!inherits(core_polygon_df, "data.frame")) {
      stop("`core_polygon_df` must be a data.frame or NULL.", call. = FALSE)
    }
    if (!all(c("x", "y") %in% names(core_polygon_df))) {
      stop("`core_polygon_df` must contain columns `x` and `y`.", call. = FALSE)
    }
    if (!is.numeric(core_polygon_df$x) || !is.numeric(core_polygon_df$y)) {
      stop("Columns `x` and `y` in `core_polygon_df` must be numeric.", call. = FALSE)
    }
  }
  
  # logical flags
  logical_args <- list(
    use_rect_zone = use_rect_zone,
    fill_around   = fill_around,
    verbose       = verbose
  )
  
  for (nm in names(logical_args)) {
    if (!is.logical(logical_args[[nm]]) || length(logical_args[[nm]]) != 1) {
      stop(sprintf("`%s` must be a single logical value.", nm), call. = FALSE)
    }
  }
  
  
  # To wrap warnings and messages inside the function
  my_warning <- function(msg) if(verbose) warning(msg, call. = FALSE)

  
  # DEFINE THE CORE INVENTORY ZONES ----
  
  ## Check if the core polygone is a polygone (at least 3 tops) ----
  if (!is.null(core_polygon_df) && nrow(core_polygon_df) < 3) {
    my_warning("The core polygone has less than 3 tops: we did not considered it")
    core_polygon_df <- NULL
  }
  
  
  ## Define a polygon that encompasses all the trees ----
  # Either the polygon is initially given, or created from concaveman if it was NULL
  
  coords_sf <- st_as_sf(
    dplyr::bind_rows(
      trees_inv[,c("x", "y")],
      sensors[,c("x", "y")]
    ), 
    coords = c("x", "y")
  ) 
  
  if (is.null(core_polygon_df)) {
    # Create concave hull from tree and sensor points
    core_polygon_sf <- concaveman::concaveman(coords_sf, concavity = 10)[[1]]
    
  } else {
    # Create polygon from user-supplied data.frame
    core_polygon_sf <- sfheaders::sf_polygon(core_polygon_df)
  }
  
  
  ## Ensure the polygon is valid ----
  if (!st_is_valid(core_polygon_sf)) {
    
    # Try to fix invalid geometry
    core_polygon_sf <- st_make_valid(core_polygon_sf)
    
    # If the result is a GEOMETRYCOLLECTION, extract POLYGONs only
    if (any(grepl("GEOMETRYCOLLECTION", class(core_polygon_sf)))) {
      
      # Try to extract a POLYGON
      extracted <- tryCatch({
        st_collection_extract(core_polygon_sf, "POLYGON")
      }, error = function(e) {
        NULL
      })
      
      # If extraction failed or result is empty, throw error
      if (is.null(extracted) || length(extracted) == 0) {
        reason <- st_is_valid(core_polygon_sf, reason = TRUE)
        stop(paste("Could not extract a valid POLYGON from GEOMETRYCOLLECTION. Reason:", reason))
      }
      
      core_polygon_sf <- extracted
    }
    
    # Recheck that the result is valid
    if (!st_is_valid(core_polygon_sf)) {
      reason <- st_is_valid(core_polygon_sf, reason = TRUE)
      stop(paste("Polygon is still invalid after attempting to fix. Reason:", reason))
    }
    
    if (is.null(core_polygon_df)) {
      my_warning("We modify the core polygon to make it valid: check anyway you core polygon")
    }
  }
  
  
  ## Ensure that all trees and sensors are in the core polygon ----
  # use st_intersect and not st_within to also consider points in the edges of the polygon
  if (!all(st_intersects(coords_sf, core_polygon_sf, sparse = FALSE))) {
    
    ## Buffer to include edge cases ----
    # i.e. points not included because of rounding errors 
    # or maybe concaveman algorithm that missed some points
    
    # Increase progressively the buffer zone with a maximum of 1m
    buffer_dist <- 10^(-8:0)
    buffered_worked <- FALSE
    
    for (dist in buffer_dist) {
      
      # Buffer the polygon
      core_polygon_buffered_sf <- st_buffer(core_polygon_sf, dist = dist)
      
      # Check if the buffer worked
      if (all(st_intersects(coords_sf, core_polygon_buffered_sf, sparse = FALSE))) {
        buffered_worked <- TRUE
        
        # Precise with a warning the buffer
        my_warning(paste0("We added a ", dist, "m buffer around the core polygon because of difficulties of including some points..."))
        
        # Keep the buffered polygon
        core_polygon_sf <- core_polygon_buffered_sf
        
        # Stop the loop
        break
      }
    }
    
    # If it does not work even after a 1m buffer
    if (!buffered_worked) {
      stop("Some trees or sensors are outside the core polygon (even after buffering of about 1 meter). Check your core polygon data.frame or create one yourself to ensure including all the tree and sensor points (automatic algorithm failed)...")
    }
    
  }
  
  
  ## Get the final core polygon df ----
  core_polygon_df <- st_coordinates(core_polygon_sf) %>% 
    as.data.frame() %>% 
    dplyr::select(x = X, y = Y)
  
  
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
  n_cells_x <- ceiling(inv_size_x / cell_size)
  n_cells_y <- ceiling(inv_size_y / cell_size)
  
  
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
  trees_shifted <- trees_inv %>% 
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
  )
  
  
  
  # FILL THE TABLE WITH LAST VARIABLES ----
  
  ## Compute tree z ----
  trees_filled <- trees_filled %>% 
    dplyr::mutate(
      z = get_z(x, y,
                deg2rad(slope),
                deg2rad(get_bottom_azimut(aspect, north2x))),
      .after = y
    )
  
  ## Compute tree hmax ----
  trees_filled <- trees_filled %>%
    dplyr::mutate(
      hmax_m = dplyr::case_when(
        
        # Paraboloid symmetric and assymetric crowns
        # Crown maximum radius fixed at crown base height
        crown_type %in% c("P", "4P") ~ hbase_m,
        
        # Symmetric and assymetric ellipsoidal crowns 
        # with fixed mid-height maximum radius
        crown_type %in% c("E", "4E") ~ hbase_m + 0.5 * (h_m - hbase_m),
        
        # Symmetric and assymetric ellipsoidal crowns 
        # with user-defined hmax
        crown_type %in% c("2E", "8E") ~ hmax_m,
        
        # Safety fallback (should never occur if crown_type was checked)
        TRUE ~ NA_real_
      )
    )
  
  ## Add dbh column to NA if it does not exist ----
  if (!"dbh_cm" %in% names(trees_filled)) {
    trees_filled <- trees_filled %>% 
      dplyr::mutate(dbh_cm = NA_real_, .after = z)
  }
  
  ## Add crown openness / LAD columns if missing ----
  if (!"crown_openness" %in% names(trees_filled)) {
    trees_filled <- trees_filled %>% 
      dplyr::mutate(crown_openness = NA_real_)
  }
  
  if (!"crown_lad" %in% names(trees_filled)) {
    trees_filled <- trees_filled %>% 
      dplyr::mutate(crown_lad = NA_real_)
  }
  
  ## Convert as data.frame ----
  trees_filled <- trees_filled %>% as.data.frame()
  
  
  # CREATE THE GRID OF CELLS ----
  cells <- tidyr::expand_grid(
    c = 0:(n_cells_x - 1),
    r = 0:(n_cells_y - 1)
  ) %>%
    dplyr::mutate(
      id_cell = n_cells_x * r + c + 1,
      
      x_center = cell_size * (c + 0.5),
      y_center = cell_size * (n_cells_y - r - 0.5),
      
      z_center = get_z(
        x_center,
        y_center,
        deg2rad(slope),
        deg2rad(get_bottom_azimut(aspect, north2x))
      )
    ) %>%
    dplyr::select(
      id_cell,
      x_center,
      y_center,
      z_center,
      row = r,
      col = c
    ) %>%
    dplyr::arrange(id_cell) %>%
    as.data.frame()
  
  
  # SHIFT THE SENSOR DATAFRAME
  if (!is.null(sensors)) {
    sensors$x <- sensors$x + shift_x
    sensors$y <- sensors$y + shift_y
    
    sensors$z <- get_z(
      sensors$x, sensors$y,
      deg2rad(slope), 
      deg2rad(get_bottom_azimut(aspect, north2x))
    ) + sensors$h_m
  }
  
  
  # OUTPUT A S3 OBJECT  ----

  ## Create the object ----
  stand <- list(
    "trees" = trees_filled,
    "sensors" = sensors,
    "cells" = cells,
    "core_polygon" = list(
      "df" = inv_zone_shifted_df,
      "sf" = inv_zone_shifted_sf,
      "use_rect_zone" = use_rect_zone
    ),
    "transform" = list(
      "core_area_ha" = core_area_ha, 
      "core_batot_m2ha" = batot_target_m2ha,
      "fill_around" = fill_around,
      "n_added_tree" = i,
      "new_area_ha" = plot_area_ha,
      "new_batot_m2ha" = current_batot_m2ha,
      "shift_x" = shift_x,
      "shift_y" = shift_y
    ),
    "geometry" = list("cell_size" = cell_size,
                      "n_cells_x" = n_cells_x, 
                      "n_cells_y" = n_cells_y,
                      "slope" = slope,
                      "aspect" = aspect,
                      "north2x" = north2x
    )
  )
  
  class(stand) <- c("sl_stand", "list")
  
  ## Validate the object ----
  validate_sl_stand(stand)
  if (verbose) message("SamsaRaLight stand successfully created.")
  
  return(stand)
}


#' Validate a SamsaRaLight stand object
#'
#' This function checks the internal consistency and structure of an object of
#' class \code{"sl_stand"}, as returned by \link{create_sl_stand}. It verifies that
#' all required components are present and correctly formatted, that the
#' embedded tree inventory conforms to the rules enforced by
#' \link{check_inventory}, and that all trees and sensors are within the stand limits.
#'
#' @param x An object expected to be of class \code{"sl_stand"}.
#'
#' @details
#' The following validations are performed:
#' \itemize{
#'   \item The object inherits from class \code{"sl_stand"}.
#'   \item The top-level components \code{trees}, \code{sensors}, \code{cells},
#'   \code{core_polygon}, \code{transform}, and \code{geometry} are present.
#'   \item The \code{trees} data.frame passes \link{check_inventory}.
#'   \item The \code{sensors} data.frame passes \link{check_sensors}.
#'   \item The \code{cells} data.frame contains columns \code{x_center},
#'   \code{y_center}, \code{z_center}, and \code{id_cell}.
#'   \item The \code{geometry} list contains \code{cell_size}, \code{n_cells_x},
#'   \code{n_cells_y}, \code{slope}, \code{aspect}, and \code{north2x}.
#'   \item All trees and sensors lie within the bounds of the rectangular stand.
#' }
#'
#' @return Invisibly returns \code{TRUE} if the stand is valid.
#'
#' @seealso
#' \link{create_sl_stand}, \link{check_inventory}, \link{check_sensors}, \link{run_sl}
#'
#' @examples
#' \dontrun{
#' data_prenovel <- SamsaRaLight::data_prenovel
#' stand <- create_sl_stand(data_prenovel$trees, cell_size = 5,
#'                          slope = 10, aspect = 180, north2x = 0)
#'
#' validate_sl_stand(stand)  # returns TRUE invisibly
#' }
#'
#' @keywords internal
validate_sl_stand <- function(x) {
  
  # Class check
  if (!inherits(x, "sl_stand")) {
    stop("Object is not a `sl_stand`.", call. = FALSE)
  }
  
  # Top-level components
  required <- c("trees", "sensors", "cells", "core_polygon", "transform", "geometry")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("sl_stand is missing element(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  # Trees & sensors format
  check_inventory(x$trees, verbose = FALSE)
  check_sensors(x$sensors, verbose = FALSE)
  
  # Cells
  if (!is.data.frame(x$cells)) stop("`cells` must be a data.frame.", call. = FALSE)
  required_cells <- c("x_center","y_center","z_center","id_cell")
  if (!all(required_cells %in% names(x$cells))) {
    stop("`cells` is malformed. Missing column(s): ",
         paste(setdiff(required_cells, names(x$cells)), collapse = ", "), call. = FALSE)
  }
  
  # Geometry
  geom <- x$geometry
  if (!is.list(geom)) stop("`geometry` must be a list.", call. = FALSE)
  required_geom <- c("cell_size","n_cells_x","n_cells_y","slope","aspect","north2x")
  if (!all(required_geom %in% names(geom))) {
    stop("`geometry` is incomplete. Missing element(s): ",
         paste(setdiff(required_geom, names(geom)), collapse = ", "), call. = FALSE)
  }
  
  # Stand limits
  x_max <- geom$cell_size * geom$n_cells_x
  y_max <- geom$cell_size * geom$n_cells_y
  if (any(x$trees$x < 0 | x$trees$x > x_max |
          x$trees$y < 0 | x$trees$y > y_max)) {
    stop("Some trees are outside the stand limits [0,", x_max, "] x [0,", y_max, "].",
         call. = FALSE)
  }
  if (!is.null(x$sensors) && nrow(x$sensors) > 0) {
    if (any(x$sensors$x < 0 | x$sensors$x > x_max |
            x$sensors$y < 0 | x$sensors$y > y_max)) {
      stop("Some sensors are outside the stand limits [0,", x_max, "] x [0,", y_max, "].",
           call. = FALSE)
    }
  }
  
  invisible(TRUE)
}