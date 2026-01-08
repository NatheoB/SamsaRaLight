# Function to plot SamsaRaLight object output

@param plot.trees Character string indicating the filling variables of
the trees. Set it to NULL if you do not want to plot the trees

## Usage

``` r
plot_sl_output(
  sl_output,
  trees.only_inventoried = FALSE,
  trees.border.species = FALSE,
  trees.fill = "e",
  trees.fill.palette = c("viridis", "light", "base"),
  trees.fill.inverse = FALSE,
  trees.fill.limits = NULL,
  cells.border = FALSE,
  cells.fill = "e",
  cells.fill.palette = c("light", "viridis", "base"),
  cells.fill.inverse = FALSE,
  cells.fill.limits = NULL,
  sensors.plot = FALSE
)
```
