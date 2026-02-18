# Plot a compass showing stand orientation and slope aspect

This function creates a small graphical compass illustrating the
orientation of the stand coordinate system relative to North, together
with the slope aspect direction if a slope is defined.

## Usage

``` r
plot_orientation_compass(north2x, slope, aspect, size = 1)
```

## Arguments

- north2x:

  Numeric. Clockwise angle (in degrees) from geographic North to the
  X-axis of the planar coordinate system.

- slope:

  Numeric. Slope angle in degrees. If `NA` or equal to 0, the aspect
  arrow is not displayed.

- aspect:

  Numeric. Aspect of the slope in degrees, measured clockwise from
  geographic North. Only used if `slope > 0`.

- size:

  Numeric. Scaling factor controlling the size of the compass arrows
  (default = 1).

## Value

A `ggplot2` object representing the orientation compass.

## Details

Cardinal directions (N, E, S, W) are drawn according to the `north2x`
parameter, which defines the clockwise rotation (in degrees) from
geographic North to the X-axis of the planar coordinate system. The
North direction is highlighted in red.

When `slope > 0`, the slope aspect (i.e. the downslope direction, given
as a clockwise angle from North) is shown as a dashed blue arrow.
