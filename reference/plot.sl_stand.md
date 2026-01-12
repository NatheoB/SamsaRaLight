# Plot a SamsaRaLight virtual stand

This function plots a virtual forest stand (`sl_stand`) produced by
SamsaRaLight. It can display a top-down view with tree crowns or a
side/top view with cells and trees.

## Usage

``` r
# S3 method for class 'sl_stand'
plot(x, ..., top_down = FALSE, only_inv = FALSE, add_sensors = TRUE)
```

## Arguments

- x:

  An object of class `sl_stand`.

- ...:

  Additional arguments passed to lower-level plotting functions.

- top_down:

  Logical, if TRUE, creates a top-down view with multiple directions
  (south, north, west, east).

- only_inv:

  Logical, if TRUE, plot only trees from the initial inventory (i.e. not
  trees added to fill around the core polygon)

- add_sensors:

  Logical; if TRUE (default), sensors are drawn on the plot. In top-down
  mode, sensors are shown as segment from ground to their height; in map
  view, sensors are drawn as squares on the ground.

## Value

A `ggplot` object representing the stand.

## Details

For the sake of the representation in top-down plot, z are offset such
as minimum altitude tree is at Y-axis height = 0
