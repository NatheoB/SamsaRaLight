# Axis-align an inventory polygon and rotate coordinates accordingly

This internal function ensures that the inventory zone is an
axis-aligned rectangle. If the input polygon is not already a rectangle,
the minimum-area enclosing rectangle is computed and rotated to align
with the X and Y axes using the minimum possible rotation.

## Usage

``` r
create_aarect_inventory(core_polygon_df, trees_inv, north2x, sensors = NULL)
```

## Arguments

- core_polygon_df:

  A validated data.frame defining the inventory polygon with columns `x`
  and `y`

- trees_inv:

  A data.frame containing tree coordinates (`x`, `y`)

- north2x:

  Numeric. Clockwise angle from North to the X-axis (degrees)

- sensors:

  Optional data.frame containing sensor coordinates (`x`, `y`), or
  `NULL`

## Value

A list with elements:

- `core_polygon_df`:

  Axis-aligned rectangular polygon

- `trees_inv`:

  Updated tree data.frame

- `sensors`:

  Updated sensor data.frame (or NULL)

- `north2x`:

  Updated orientation angle

- `rotation_rad`:

  Applied rotation angle (radians)

## Details

Tree and sensor coordinates are rotated consistently.
