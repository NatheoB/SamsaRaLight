#' Plot a from-above view of a tree inventory
#' 
#' Visualizes a forest inventory as ellipses representing tree crowns.
#' 
#' @param trees_inv A data.frame of trees that passed \link{check_inventory}.
#' @param core_polygon_df Optional data.frame defining the core inventory polygon.
#'   Must contain columns \code{x} and \code{y}.
#' @param show_id Logical; if TRUE (default), displays tree identifiers at crown centers.
#' 
#' @details Because the north2x variable is unknown, trees are plotted as circles 
#' by considering the mean radius on the four cardinals.
#' 
#' @return A ggplot object displaying the trees in a from-above view.
#' 
#' @importFrom ggplot2 ggplot aes geom_text coord_equal theme_minimal theme xlab ylab
#' @importFrom ggforce geom_circle
#' 
#' @export
plot_inventory <- function(trees_inv, core_polygon_df = NULL, show_id = TRUE) {
  
  # Ensure coordinates are correct
  needs_conversion_trees <- check_coordinates(trees_inv, verbose = F)
  
  if (needs_conversion_trees) {
    stop(
      "Geographic coordinates (`lon`, `lat`) detected in `trees`: ",
      "convert the coordinates using `create_xy_from_lonlat()`."
    )
  }
  
  if (!is.null(core_polygon_df)) {
    needs_conversion_polygon <- check_coordinates(core_polygon_df, verbose = F)
    
    if (needs_conversion_polygon) {
      stop(
        "Geographic coordinates (`lon`, `lat`) detected in `core_polygon_df`: ",
        "convert the coordinates using `create_xy_from_lonlat()`."
      )
    }
  }
  
  # Plot the trees with the inventory zone if specified
  plt <- ggplot()
    
  # CORE POLYGON
  if (!is.null(core_polygon_df)) {  
    plt <- plt +
      geom_polygon(data = core_polygon_df,
                   mapping = aes(x = x, y = y),
                   fill = "yellow", color = "black", alpha = 0.5)
  }
  
  # TREES 
  plt <- plt +
    geom_circle(
      data = trees_inv,
      mapping = aes(
        x0 = x,
        y0 = y,
        r = (re_m + rw_m + rn_m + rs_m) / 4,
        fill = species
      ),
      color = "black",
      alpha = 0.6
    )
  
  
  # Add tree id labels if requested
  if (show_id) {
    plt <- plt +
      geom_text(
        data = trees_inv,
        aes(x = x, y = y, label = id_tree),
        size = 2,
        color = "black"
      )
  }
  
  plt +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    xlab("") + ylab("")
  
}