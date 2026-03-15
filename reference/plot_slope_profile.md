# Plot a slope profile indicator

This function draws a small schematic profile representing terrain slope
as an inclined line relative to a horizontal reference.

## Usage

``` r
plot_slope_profile(slope, size = 1)
```

## Arguments

- slope:

  Numeric. Slope angle in degrees. If `NA` or `0`, only the horizontal
  reference is shown.

- size:

  Numeric. Scaling factor controlling the size of the graphic (default =
  1).

## Value

A `ggplot2` object representing the slope profile.

## Details

It is intended to be used as a companion graphic to stand orientation
plots (e.g. compass), for visualizing slope magnitude.
