#' Plot a from-above view of a tree inventory
#' 
#' Visualizes a forest inventory as ellipses representing tree crowns.
#' 
#' @param trees_inv A data.frame of trees that passed \link{check_inventory}.
#' @param show_id Logical; if TRUE (default), displays tree identifiers at crown centers.
#' 
#' @return A ggplot object displaying the trees as ellipses in a from-above view.
#' 
#' @importFrom ggplot2 ggplot aes geom_text coord_equal theme_minimal theme xlab ylab
#' @importFrom ggforce geom_ellipse
#' 
#' @export
plot_inventory <- function(trees_inv, show_id = TRUE) {
  
  # Ensure inventory is correct
  check_inventory(trees_inv, verbose = F)
  
  
  # Plot the trees
  plt <- ggplot() +
    
    geom_ellipse(
      data = trees_inv,
      mapping = aes(
        x0 = x,
        y0 = y,
        a = (re_m + rw_m) / 2,
        b = (rn_m + rs_m) / 2,
        angle = 0,
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