# Make the inventory polygon rectangle and optionally axis-align it

The minimum-area enclosing rectangle is computed, and if
`rotate_axisaligned = TRUE`, align the rectangle with the X and Y axes,
and consequently, coordinates of trees (and sensors) are rotated.

## Usage

``` r
create_rect_inventory(
  core_polygon_df,
  trees,
  north2x,
  sensors = NULL,
  rotate_axisaligned = TRUE
)
```

## Arguments

- core_polygon_df:

  A validated data.frame defining the inventory polygon with columns `x`
  and `y`.

- trees:

  A data.frame containing tree coordinates (`x`, `y`).

- north2x:

  Numeric. Clockwise angle from North to the X-axis (degrees).

- sensors:

  Optional data.frame containing sensor coordinates (`x`, `y`), or
  `NULL`.

- rotate_axisaligned:

  Logical. If `TRUE`, rotate coordinates to align rectangle with axes.

## Value

A list with elements:

- `core_polygon_df`:

  Rectangular polygon (rotated axis-aligned if requested).

- `trees`:

  Updated tree data.frame (rotated if requested).

- `sensors`:

  Updated sensor data.frame (rotated if requested, or NULL).

- `north2x`:

  Updated orientation angle (degrees).

- `rotation`:

  Applied rotation angle in degrees. 0 if no rotation.

## Details

Tree and sensor coordinates are rotated consistently.
