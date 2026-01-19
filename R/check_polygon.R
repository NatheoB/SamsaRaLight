#' Check and validate a polygon defined by vertices
#'
#' This function converts a data.frame of polygon vertices into an sf POLYGON
#' and checks its validity. If the polygon is invalid, it attempts to fix it.
#'
#' @param polygon_df A data.frame with columns x and y defining polygon vertices
#' @param trees_inv A data.frame with one row per tree.
#'   See \link{check_inventory} for the required structure and validated columns.
#' @param sensors Optional data.frame defining position and height of the sensor within the stand.
#'   See \link{check_sensors} for the required structure and validated columns.
#' @param verbose Logical. If TRUE, warnings are printed
#'
#' @return A data.frame of polygon vertices (x, y):
#'   - unchanged if valid
#'   - modified if fixed (with a warning)
#'
#' @export
check_polygon <- function(polygon_df, trees_inv, sensors = NULL, verbose = TRUE) {
  
  # Check trees_inv format ----
  if (!check_inventory(trees_inv, verbose = FALSE)) {
    stop("`trees_inv` must be a data.frame verified by check_inventory().", call. = FALSE)
  }
  
  # Check sensors format ----
  if (! (is.null(sensors) || check_sensors(trees_inv, verbose = FALSE)) ) {
    stop("`sensors` must be NULL or a data.frame verified by check_sensors().", call. = FALSE)
  }
  
  # Check data.frame format of polygon ----
  if (!inherits(polygon_df, "data.frame")) {
    stop("`polygon_df` must be a data.frame.", call. = FALSE)
  }
  
  if (!all(c("x", "y") %in% names(polygon_df))) {
    stop("`polygon_df` must contain columns `x` and `y`.", call. = FALSE)
  }
  
  if (!is.numeric(polygon_df$x) || !is.numeric(polygon_df$y)) {
    stop("Columns `x` and `y` in `polygon_df` must be numeric.", call. = FALSE)
  }
  
  
  # Check if the core polygone is a polygone (at least 3 tops) ----
  if (nrow(polygon_df) < 3) {
    stop("The polygon has less than 3 vertices and cannot be formed.", call. = FALSE)
  }
  
  
  # Create polygon from user-supplied data.frame ----
  core_polygon_sf <- sfheaders::sf_polygon(polygon_df)
  
  
  # Ensure the polygon is valid ----
  if (!sf::st_is_valid(core_polygon_sf)) {
    
    # Try to fix invalid geometry
    core_polygon_sf <- sf::st_make_valid(core_polygon_sf)
    
    # If the result is a GEOMETRYCOLLECTION, extract POLYGONs only
    if (any(grepl("GEOMETRYCOLLECTION", class(core_polygon_sf)))) {
      
      # Try to extract a POLYGON
      extracted <- tryCatch({
        sf::st_collection_extract(core_polygon_sf, "POLYGON")
      }, error = function(e) {
        NULL
      })
      
      # If extraction failed or result is empty, throw error
      if (is.null(extracted) || length(extracted) == 0) {
        reason <- sf::st_is_valid(core_polygon_sf, reason = TRUE)
        stop(
          paste("Could not extract a valid POLYGON from GEOMETRYCOLLECTION. Reason:", reason),
          call. = FALSE
        )
      }
      
      core_polygon_sf <- extracted
    }
    
    # Recheck that the result is valid
    if (!sf::st_is_valid(core_polygon_sf)) {
      reason <- sf::st_is_valid(core_polygon_sf, reason = TRUE)
      stop(
        paste("Polygon is still invalid after attempting to fix. Reason:", reason),
        call. = FALSE
      )
    }
    
    if (verbose) warning("The polygon was invalid and has been modified to make it valid.", call. = FALSE)
  }
  

  # Ensure that all trees and sensors are in the core polygon ----
  coords_sf <- st_as_sf(
    dplyr::bind_rows(
      trees_inv[,c("x", "y")],
      sensors[,c("x", "y")]
    ), 
    coords = c("x", "y")
  ) 
  
  # use st_intersect and not st_within to also consider points in the edges of the polygon
  if (!all(st_intersects(coords_sf, core_polygon_sf, sparse = FALSE))) {
    
    # Increase progressively the buffer zone with a maximum of ~1‰ of polygon size
    # (scale-aware: works for metric and normalized coordinates)
    bbox <- st_bbox(core_polygon_sf)
    poly_scale <- sqrt((bbox$xmax - bbox$xmin)^2 +
                         (bbox$ymax - bbox$ymin)^2)
    
    # Relative buffer distances (from ~1e-8 to ~1e-3 of polygon size)
    buffer_dist <- poly_scale * 10^(-8:-3)
    buffered_worked <- FALSE
    
    for (dist in buffer_dist) {
      
      # Buffer the polygon
      core_polygon_buffered_sf <- st_buffer(core_polygon_sf, dist = dist)
      
      # Check if the buffer worked
      if (all(st_intersects(coords_sf, core_polygon_buffered_sf, sparse = FALSE))) {
        buffered_worked <- TRUE
        
        # Precise with a warning the buffer
        if (verbose)
          warning(paste0(
            "We added a buffer of ",
            signif(dist, 3),
            " (≈ ",
            signif(dist / poly_scale * 100, 3),
            "% of polygon size) around the core polygon to include edge points."
          ))
        
        # Keep the buffered polygon
        core_polygon_sf <- core_polygon_buffered_sf
        
        # Stop the loop
        break
      }
    }
    
    # If it does not work even after the maximum relative buffer
    if (!buffered_worked) {
      stop(
        "Some trees or sensors are outside the core polygon even after ",
        "relative buffering based on polygon scale. ",
        "Check your core polygon data.frame or create one yourself ",
        "to ensure including all the tree and sensor points ",
        "(automatic algorithm failed)..."
      )
    }
    
  }
  
  # Get the final polygon df ----
  polygon_checked_df <- sf::st_coordinates(core_polygon_sf) %>% 
    as.data.frame() %>% 
    dplyr::select(x = X, y = Y)
  
  
  ## Success ----
  if (verbose) message("Polygon successfully validated.")
  
  return(polygon_checked_df)
}
