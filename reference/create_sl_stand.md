# Create a rectangular virtual stand from a tree inventory

This function builds a rectangular (or square) virtual forest stand from
a user-provided tree inventory. Trees are spatially shifted so that the
inventory zone is centered within a regular grid of square cells.
Optionally, additional trees can be added around the core inventory area
to match its basal area per hectare.

## Usage

``` r
create_sl_stand(
  trees_inv,
  cell_size,
  latitude,
  slope,
  aspect,
  north2x,
  sensors = NULL,
  core_polygon_df = NULL,
  aarect_zone = FALSE,
  fill_around = FALSE,
  verbose = TRUE
)
```

## Arguments

- trees_inv:

  A data.frame with one row per tree. See
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md)
  for the required structure and validated columns.

- cell_size:

  Numeric. Side length of square cells composing the stand (meters).

- latitude:

  Numeric, latitude of the stand (degrees)

- slope:

  Numeric. Slope of the plot (degrees).

- aspect:

  Numeric. Aspect of the slope, defined as the azimuth of the downslope
  direction, clockwise from North (degrees). North = 0, East = 90, South
  = 180, West = 270.

- north2x:

  Numeric. Clockwise angle from North to the X-axis (degrees). A value
  of 0 corresponds to a Y-axis oriented toward North.

- sensors:

  Optional data.frame defining position and height of the sensor within
  the stand. See
  [check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md)
  for the required structure and validated columns.

- core_polygon_df:

  Optional data.frame defining the core inventory polygon. Must contain
  columns `x` and `y`. If `NULL`, a concave hull is automatically
  computed from tree positions.

- aarect_zone:

  Logical. If `TRUE`, the inventory zone is defined by the minimum-area
  enclosing rectangle of the core polygon with minimum rotation to
  obtain an axis-aligned rectangle inventory zone. If `FALSE`, the core
  polygon itself is used.

- fill_around:

  Logical. If `TRUE`, trees are added outside the core polygon until the
  basal area per hectare of the full stand matches that of the core
  inventory.

- verbose:

  Logical. If `TRUE` (default), messages and warnings are printed during
  processing. If `FALSE`, output is silent.

## Value

A named list with the following elements:

- `trees`:

  Data.frame of the final tree inventory, including added trees if
  `fill_around = TRUE`. The structure matches the inventory format
  validated by
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
  with additional derived variables required for ray tracing. A logical
  column `added_to_fill` indicates whether each tree originates from the
  initial inventory or was added to fill around the inventory zone.

- `core_polygon`:

  List describing the inventory zone:

  - `df`: data.frame of polygon vertices

  - `sf`: corresponding `sf` POLYGON

  - `aarect_zone`: did we used an axis-aligned rectangle inventory zone
    ?

- `transform`:

  List of transformation and filling information, including core area,
  target and final basal area, number of added trees, and applied
  spatial shifts.

- `geometry`:

  List describing stand geometry and terrain parameters (cell size,
  number of cells, slope, aspect, and orientation).

## Details

The function supports sloping terrain and coordinate system rotation,
and returns a fully prepared stand ready for use in the ray-tracing
pipeline (see
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md)).

The returned `trees` data.frame conforms to the inventory format checked
by
[check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
with the following controlled modifications:

- Tree vertical position `z` is computed from terrain slope and aspect.

- Crown maximum radius height `hmax_m` is computed when fixed by crown
  geometry conventions:

  - `"P"` and `"4P"`: `hmax_m = hbase_m`

  - `"E"` and `"4E"`: `hmax_m = hbase_m + 0.5 * (h_m - hbase_m)`

  For crown types `"2E"` and `"8E"`, `hmax_m` must be provided by the
  user.

- If missing, column `dbh_cm` is added and filled with `NA`.

- If missing, crown interception properties (e.g. `crown_openness`,
  `crown_lad`) are added using default values.

The function ensures that all trees fall within the core inventory
polygon, applying small buffers if necessary to handle numerical
precision issues. Invalid polygons are automatically repaired when
possible.

When `fill_around = TRUE`, trees are randomly sampled from the original
inventory and positioned outside the core polygon until the target basal
area per hectare is reached for the full rectangular stand.

## Examples

``` r
if (FALSE) { # \dontrun{
data_prenovel <- SamsaRaLight::data_prenovel
trees_inv <- data_prenovel$trees

stand <- create_sl_stand(
  trees_inv = trees_inv,
  cell_size = 5,
  latitude = 46,
  slope = 10,
  aspect = 180,
  north2x = 0,
  aarect_zone = TRUE,
  fill_around = FALSE,
  verbose = TRUE
)

head(stand$trees)
} # }
```
