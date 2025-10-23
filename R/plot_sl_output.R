#' @import ggplot2 ggforce
#'
#' @export
#'
plot_sl_output <- function(sl_output, 
                           plot.trees = TRUE) {
  
  cell_size <- sl_output$input$info$cell_size
  n_cells_x <- sl_output$input$info$n_cells_x
  n_cells_y <- sl_output$input$info$n_cells_y
  
  plt <- sl_output$output$cells %>% 
    dplyr::mutate(
      pacl_group = dplyr::case_when(
        pacl > 0.5 ~ "> 50%",
        pacl > 0.25 ~ "25 - 50%",
        pacl > 0.125 ~ "12.5 - 25%",
        pacl > 0.0625 ~ "6.25 - 12.5%",
        TRUE ~ "< 6.25%"
      ),
      pacl_group = factor(pacl_group,
                          levels = c("< 6.25%", "6.25 - 12.5%",
                                     "12.5 - 25%", "25 - 50%",
                                     "> 50%"))) %>% 
    
    ggplot() +
    
    geom_tile(aes(x = x_center, y = y_center,
                  fill = pacl_group), color = "black") +
    scale_fill_grey(start = 0.2, end = 1,
                    drop = FALSE) +
    labs(fill = "Cells enlightment") +
    
    
    coord_equal() +
    scale_x_continuous(breaks = seq(0, n_cells_x*cell_size,
                                    by = cell_size),
                       labels = round(seq(0, n_cells_x*cell_size,
                                          by = cell_size),
                                      digits = 1)) +
    scale_y_continuous(breaks = seq(0, n_cells_y*cell_size,
                                    by = cell_size),
                       labels = round(seq(0, n_cells_y*cell_size,
                                          by = cell_size),
                                      digits = 1)) +
    
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right") +
    xlab("") + ylab("")
  
  # Plot trees if precised  
  if (plot.trees) {
    
    plt <- plt + 
      ggnewscale::new_scale_fill() +
      geom_ellipse(data = sl_output$input$trees, 
                   mapping = aes(x0 = x, y0 = y, 
                                 a = pmax(re_m, rw_m),
                                 b = pmax(rn_m, rs_m),
                                 angle = 0,
                                 fill = species)) +
      labs(fill = "Species")
  }
  
  plt
}