# Check the format and validity of a tree inventory data.frame

This function checks whether a tree inventory data.frame is correctly
formatted to be used as input for the ray-tracing model. It verifies the
presence, type, and validity of mandatory and optional variables
describing the geometry and attributes of trees within a forest stand.

## Usage

``` r
check_inventory(trees_inv, verbose = TRUE)
```

## Arguments

- trees_inv:

  A data.frame with one row per tree and the following columns:

  id_tree

  :   Unique identifier of the tree (numeric or character, no
      duplicates).

  x

  :   X position of the tree within the stand (numeric, meters, planar
      coordinates). Optional if `lon` and `lat` are provided.

  y

  :   Y position of the tree within the stand (numeric, meters, planar
      coordinates). Optional if `lon` and `lat` are provided.

  lon

  :   Longitude of the tree location (numeric, decimal degrees).
      Optional; if provided with `lat`, coordinates are converted
      internally to a UTM planar coordinate system in
      [`create_sl_stand()`](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md).

  lat

  :   Latitude of the tree location (numeric, decimal degrees).
      Optional; see `lon`.

  species

  :   Species name (character).

  dbh_cm

  :   Diameter at breast height (1.3 m, in cm).

  crown_type

  :   Type of crown geometry. One of `"E"`, `"P"`, `"2E"`, `"8E"`, or
      `"4P"`.

  h_m

  :   Total tree height (numeric, meters).

  hbase_m

  :   Height of the crown base (numeric, meters).

  hmax_m

  :   Height of the maximum crown radius (numeric, meters). Required
      only if at least one tree has a crown type `"2E"` or `"8E"`. For
      other crown types, the column is optional and the value is
      internally computed.

  rn_m

  :   Crown radius toward North (numeric, meters).

  rs_m

  :   Crown radius toward South (numeric, meters).

  re_m

  :   Crown radius toward East (numeric, meters).

  rw_m

  :   Crown radius toward West (numeric, meters).

  crown_openness

  :   Crown openness (unitless), optional if turbid medium interception.

  crown_lad

  :   Leaf Area Density (m² m⁻³), optional if porous envelope
      interception.

- verbose:

  Logical; if `TRUE`, informative messages and warnings are printed.

## Value

Invisibly returns `TRUE` if all checks pass.

## Details

The function performs the following checks and validations:

- 1:

  Ensures `trees_inv` is a non-empty data.frame with all required
  columns.

- 2:

  Checks that `id_tree` values are unique.

- 3:

  Checks that either planar coordinates (`x`, `y`) or geographic
  coordinates (`lon`, `lat`) are provided. If only `lon` and `lat` are
  supplied, they are converted to planar UTM coordinates in
  [`create_sl_stand()`](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md).
  The
  [`plot_inventory()`](https://natheob.github.io/SamsaRaLight/reference/plot_inventory.md)
  function requires planar coordinates and cannot be used before this
  conversion.

- 4:

  Validates numeric columns (`dbh_cm`, `h_m`, `hbase_m`, `rn_m`, `rs_m`,
  `re_m`, `rw_m`) are numeric and non-negative.

- 5:

  Verifies that `crown_type` values are one of `"E"`, `"P"`, `"2E"`,
  `"8E"`, or `"4P"`.

- 6:

  Ensures crown radii are present according to crown type.

- 7:

  Checks that `hmax_m` is provided when required and lies between
  `hbase_m` and `h_m`.

- 8:

  Ensures `hbase_m < h_m`.

- 9:

  Verifies that each tree has at least one crown interception property
  defined.

- 10:

  Provides informative error messages and warnings for all invalid
  conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using the example dataset from the package
data_prenovel <- SamsaRaLight::data_prenovel
trees <- data_prenovel$trees

# Check the inventory
check_inventory(trees)

# Quiet mode
check_inventory(trees, verbose = FALSE)
} # }
```
