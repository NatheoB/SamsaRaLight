#' Plot a from-above view of a tree inventory
#' 
#' Visualizes a forest inventory as ellipses representing tree crowns.
#' Taller trees are drawn above shorter ones. Optionally colors trees by species.
#' 
#' @param trees_inv A data.frame of trees that passed \link{check_inventory}.
#'   Must contain at least:
#'   \describe{
#'     \item{\code{x}, \code{y}}{Tree coordinates (meters).}
#'     \item{\code{rn_m}, \code{rs_m}, \code{re_m}, \code{rw_m}}{Crown radii (meters).}
#'     \item{\code{h_m}}{Total height (meters) for plotting order.}
#'     \item{\code{id_tree}}{Optional if show_id is FALSE; Tree identifier used for labeling.}
#'     \item{\code{species}}{Optional; species name (character) for coloring.}
#'   }
#' @param transparency Logical; if TRUE (default), uses alpha = 0.6, otherwise alpha = 1 (opaque).
#' @param show_id Logical; if TRUE (default), displays tree identifiers at crown centers.
#' 
#' @details Time for plotting can be long because trees are plotted one by one in order to respect height 
#'  layering for both crowns and labels.
#'
#' @return A ggplot object displaying the trees as ellipses in a from-above view.
#' 
#' @importFrom ggplot2 ggplot aes geom_text coord_equal theme_minimal theme xlab ylab
#' @importFrom ggforce geom_ellipse
#' 
#' @export
plot_inventory <- function(trees_inv, transparency = TRUE, show_id = TRUE) {
  
  # Ensure inventory is correct
  check_inventory(trees_inv, verbose = F)
  
  # Set alpha based on transparency argument
  alpha_val <- ifelse(transparency, 0.6, 1)
  
  # Order trees by height (smaller first, taller on top)
  trees_plot <- trees_inv[order(trees_inv$h_m), ]
  
  
  # Plot the trees with label in the height order
  # Mandatory to plot row-wise for plotting labels also in height order
  # And because if considering species, plotting order is by heights but grouped by species (not over all trees)
  plt <- ggplot()
  
  for (i in seq_len(nrow(trees_plot))) {
    
    plt <- plt +
      geom_ellipse(
        data = trees_plot[i,],
        mapping = aes(
          x0 = x,
          y0 = y,
          a = (re_m + rw_m) / 2,
          b = (rn_m + rs_m) / 2,
          angle = 0,
          fill = species
        ),
        color = "black",
        alpha = alpha_val
      )
    
    
    # Add tree id labels if requested
    if (show_id) {
      plt <- plt +
        geom_text(
          data = trees_plot[i,],
          aes(x = x, y = y, label = id_tree),
          size = 2,
          color = "black"
        )
    }
    
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