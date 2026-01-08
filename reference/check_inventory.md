# Check the format and validity of a tree inventory data.frame

This function checks whether a tree inventory data.frame is correctly
formatted to be used as input for the ray-tracing model. It verifies the
presence, type, and validity of mandatory and optional variables
describing the geometry and attributes of trees within a forest stand.

## Usage

``` r
check_inventory(tree_inv, verbose = TRUE)
```

## Arguments

- tree_inv:

  A data.frame with one row per tree and the following columns:

  - id_treeUnique identifier of the tree (numeric or character, no
    duplicates)

  - xX position of the tree within the stand (numeric, meters)

  - yY position of the tree within the stand (numeric, meters)

  - speciesOPTIONAL – species name (character)

  - dbh_cmDiameter at breast height (1.3 m, in cm)

  - crown_typeType of crown geometry. One of `"E"`, `"P"`, `"2E"`,
    `"8E"`, or `"4P"`

  - h_mTotal tree height (numeric, meters)

  - hbase_mHeight of the crown base (numeric, meters)

  - hmax_mHeight of the maximum crown radius (numeric, meters). Required
    only if at least one tree has a crown type `"2E"` or `"8E"`. For
    other crown types, the value is ignored and internally recomputed.

  - rn_mCrown radius toward North (numeric, meters)

  - rs_mCrown radius toward South (numeric, meters)

  - re_mCrown radius toward East (numeric, meters)

  - rw_mCrown radius toward West (numeric, meters)

  - crown_opennessOPTIONAL – Crown openness (unitless), mandatory for
    porous envelope interception

  - crown_ladOPTIONAL – Leaf Area Density (m² m⁻³), mandatory for turbid
    medium interception

- verbose:

  Logical; if `TRUE`, informative messages and warnings are printed.

## Value

Invisibly returns `TRUE` if all checks pass.

## Details

The function performs the following checks and validations:

- Ensures `tree_inv` is a non-empty data.frame with all required
  columns.

- Checks that `id_tree` values are unique.

- Validates numeric columns (`x`, `y`, `dbh_cm`, `h_m`, `hbase_m`,
  `rn_m`, `rs_m`, `re_m`, `rw_m`) are numeric and non-negative.

- Verifies that `crown_type` values are one of `"E"`, `"P"`, `"2E"`,
  `"8E"`, or `"4P"`.

- Ensures crown radii are present according to crown type. For symmetric
  crowns (`"E"`, `"P"`, `"2E"`), a single effective radius is computed
  if the four directional radii differ. For asymmetric crowns (`"8E"`,
  `"4P"`), all four directional radii must be provided.

- Checks that `hmax_m` is provided when required by crown type (`"2E"`
  or `"8E"`) and lies between `hbase_m` and `h_m`.

- Ensures `hbase_m` \< `h_m`.

- Verifies that each tree has at least one crown interception property
  defined (`crown_openness` or `crown_lad`).

- Provides informative error messages and warnings for all invalid
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
