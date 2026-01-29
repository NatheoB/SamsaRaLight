#' Plot a SamsaRaLight output
#'
#' Visualize ground light and tree energy metrics for a \code{sl_output} object.
#'
#' @param x An object of class \code{sl_output}, returned by \code{run_sl()}.
#' @param ... Additional arguments passed to lower-level plotting functions.
#' @param what_trees Character; which tree metric to plot. Choices are:
#'   \describe{
#'     \item{"compet"}{Light competition index (LCI), reversed viridis scale.}
#'     \item{"intercepted"}{Intercepted energy (MJ).}
#'     \item{"potential"}{Potential intercepted energy (MJ).}
#'   }
#'   Default is "compet".
#' @param what_cells Character; which cell (ground) metric to plot. Choices are:
#'   \describe{
#'     \item{"relative"}{Proportion of above canopy light (PACL)).}
#'     \item{"absolute"}{Energy on the ground (MJ/m2).}
#'   }
#'   Default is "relative".
#' @param show_trees Logical; whether to display trees on top of the ground light map. Default is TRUE.
#' @param direct_energy Logical or NULL. 
#'   If NULL (default), total radiation outputs are plotted (direct + diffuse).
#'   If TRUE, only direct radiation components are plotted.
#'   If FALSE, only diffuse radiation components are plotted.
#'   This option requires \code{detailed_output = TRUE} when running the simulation
#' 
#' 
#' @return A ggplot object.
#' 
#' @examples
#' \dontrun{
#' plot(sl_output_object)
#' plot(sl_output_object, what_trees = "potential", what_cells = "absolute", show_trees = FALSE)
#' }
#' 
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradient scale_x_continuous
#'   scale_y_continuous coord_equal labs theme element_text element_blank guide_colorbar scale_fill_viridis_c
#' @importFrom ggplot2 xlab ylab theme_minimal
#' @importFrom scales pretty_breaks
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggforce geom_ellipse
#' 
#' @export
#' @method plot sl_output
#' 
plot.sl_output <- function(x, ...,
                           what_trees = c("compet", "intercepted", "potential"),
                           what_cells = c("relative", "absolute"),
                           show_trees = TRUE,
                           direct_energy = NULL) {
  
  # Rounding error epsilon (for 0-1 limits)
  epsilon <- 1e-4
  
  # ---- Checks ----
  stopifnot(inherits(x, "sl_output"))
  stopifnot(is.logical(show_trees), length(show_trees) == 1)
  stopifnot(is.null(direct_energy) || (is.logical(direct_energy) && length(direct_energy) == 1))
  
  what_trees <- match.arg(what_trees)
  what_cells <- match.arg(what_cells)
  
  # ---- Create base plotting tables ----
  sl_stand  <- x$input$sl_stand
  cells_out <- x$output$light$cells
  trees_out <- x$output$light$trees
  
  cells_plot <- merge(sl_stand$cells, cells_out, by = "id_cell")
  trees_plot <- merge(sl_stand$trees, trees_out, by = "id_tree")
  
  
  # ---- Energy type ----
  energy_suffix <- ""
  
  if (!is.null(direct_energy)) {
    
    # Check that detailed output exists
    if (is.null(x$input$params$detailed_output) || !x$input$params$detailed_output) {
      stop(
        "direct_energy was requested, but detailed_output = FALSE.\n",
        "Re-run the simulation with detailed_output = TRUE."
      )
    }
    
    energy_suffix <- if (direct_energy) "_direct" else "_diffuse"
  }
  
  
  # ---- Choose variables ----
  cell_var <- paste0(
    switch(
      what_cells,
      "relative" = "pacl",
      "absolute" = "e"
    ),
    energy_suffix
  )
  
  tree_var <- paste0(
    switch(
      what_trees,
      "compet"      = "lci",
      "intercepted" = "e",
      "potential"   = "epot"
    ),
    energy_suffix
  )
  
  
  # ---- Automatic legend labels ----
  tree_label <- switch(what_trees,
                       "compet" = "TREE\nLight competition index (LCI)",
                       "potential" = "TREE\nPotential intercepted energy (MJ)",
                       "intercepted" = "TREE\nIntercepted energy (MJ)")
  
  cell_label <- switch(what_cells,
                       "relative" = "CELL\nProportion of above canopy light (PACL)",
                       "absolute" = "CELL\nEnergy on the ground (MJ/m2)")
  
  
  # ---- Base: ground light ----
  plt <- ggplot() + coord_equal()
  
  if (what_cells == "relative") {
    plt <- plt +
      geom_raster(
        data = cells_plot,
        aes(x = x_center, y = y_center, fill = .data[[cell_var]])
      ) +
      scale_fill_gradient(
        name = cell_label,
        limits = c(0 - epsilon, 1 + epsilon),
        low = "black", high = "white",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          label.theme = element_text(angle = 45, hjust = 1)
        )
      )
  } else {
    plt <- plt +
      geom_raster(
        data = cells_plot,
        aes(x = x_center, y = y_center, fill = .data[[cell_var]])
      ) +
      scale_fill_gradient(
        name = cell_label,
        low = "black", high = "white",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          label.theme = element_text(angle = 45, hjust = 1)
        )
      )
  }
  
  # ---- Trees ----
  if (show_trees) {
    plt <- plt + ggnewscale::new_scale_fill()
    
    # Order trees by height
    trees_plot <- trees_plot[order(trees_plot$h_m), ]
    
    trees_plot$zval <- trees_plot[[tree_var]]
    
    
    for (i in seq_len(nrow(trees_plot))) {
      plt <- plt +
        geom_ellipse(
          data = trees_plot[i, ],
          aes(
            x0 = x,
            y0 = y,
            a = (re_m + rw_m)/2,
            b = (rn_m + rs_m)/2,
            angle = 0,
            fill = zval
          ),
          color = "black",
          linewidth = 0.3,
          alpha = 0.9
        )
    }
    
    plt <- plt +
      scale_fill_viridis_c(
        name = tree_label,
        limits = c(
          ifelse(what_trees == "compet", 0 - epsilon, NA),
          ifelse(what_trees == "compet", 1 + epsilon, NA)
        ),
        direction = ifelse(what_trees == "compet", -1, 1),
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          label.theme = element_text(angle = 45, hjust = 1)
        )
      )
  }
  
  # ---- Axes & theme ----

  energy_label <- if (is.null(direct_energy)) {
    "Total radiations"
  } else if (direct_energy) {
    "Direct radiations"
  } else {
    "Diffuse radiations"
  }
  
  xbreaks <- scales::pretty_breaks(n = 7)(seq(0, sl_stand$geometry$n_cells_x * sl_stand$geometry$cell_size, by = sl_stand$geometry$cell_size))
  ybreaks <- scales::pretty_breaks(n = 7)(seq(0, sl_stand$geometry$n_cells_y * sl_stand$geometry$cell_size, by = sl_stand$geometry$cell_size))
  
  plt <- plt +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) +
    xlab("") + ylab("") +
    labs(title = "SamsaRaLight output",
         subtitle = energy_label) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title.align = 0.5
    )
  
  plt
}
