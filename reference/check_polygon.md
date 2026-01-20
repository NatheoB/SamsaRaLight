# Check and validate a polygon defined by vertices

This function converts a data.frame of polygon vertices into an sf
POLYGON and checks its validity. If the polygon is invalid, it attempts
to fix it.

## Usage

``` r
check_polygon(polygon_df, trees_inv, sensors = NULL, verbose = TRUE)
```

## Arguments

- polygon_df:

  A data.frame with columns x and y defining polygon vertices

- trees_inv:

  A data.frame with one row per tree. See
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md)
  for the required structure and validated columns.

- sensors:

  Optional data.frame defining position and height of the sensor within
  the stand. See
  [check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md)
  for the required structure and validated columns.

- verbose:

  Logical. If TRUE, warnings are printed

## Value

A data.frame of polygon vertices (x, y):

- unchanged if valid

- modified if fixed (with a warning)
