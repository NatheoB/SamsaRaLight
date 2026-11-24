# Estimate the light on sensors and cells from a precomputed geometry

Estimates the light reaching the sensors and cells of a stand for given
crown leaf area densities (turbid medium) or crown openness (porous
envelope), from the geometry computed by
[precompute_sl_geometry](https://natheob.github.io/SamsaRaLight/reference/precompute_sl_geometry.md).
Results are the same as the `sensors` and `cells` outputs of
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
with the same settings (up to rounding errors), but much faster, since
the rays are not traced again. The light intercepted by trees cannot be
estimated from the geometry.

## Usage

``` r
estimate_pacl_from_geometry(
  sl_geometry,
  trees,
  turbid_medium = TRUE,
  extinction_coef = 0.5,
  clumping_factor = 1,
  trunk_interception = TRUE,
  detailed_output = FALSE
)
```

## Arguments

- sl_geometry:

  An object of class `"sl_geometry"` created by
  [precompute_sl_geometry](https://natheob.github.io/SamsaRaLight/reference/precompute_sl_geometry.md).

- trees:

  A data frame with one row per tree of the stand used to compute the
  geometry, with the column `id_tree`, and the column `crown_lad` (if
  `turbid_medium = TRUE`) or `crown_openness` (if
  `turbid_medium = FALSE`), as in `sl_stand$trees`. Other columns are
  ignored, and the order of the rows does not matter (trees are matched
  by `id_tree`).

- turbid_medium:

  Logical. If `TRUE`, crowns are represented as turbid media using
  `crown_lad`. If `FALSE`, crowns are represented as porous envelopes
  using `crown_openess`. Defaults to `TRUE`.

- extinction_coef:

  Numeric scalar. Leaf extinction coefficient controlling the
  attenuation of radiation by foliage. It represents the effective light
  attenuation per unit leaf area and is related to average leaf
  orientation. Higher values increase light interception. Defaults to
  `0.5`.

- clumping_factor:

  Numeric scalar controlling the aggregation of foliage within the crown
  volume. A value of `1` corresponds to homogeneous foliage
  distribution, values below `1` indicate clumped foliage, and values
  above `1` indicate more regular spacing. This parameter affects light
  interception in the turbid-medium model. Defaults to `1`.

- trunk_interception:

  Logical. If `TRUE`, account for interception by tree trunks. Defaults
  to `TRUE`.

- detailed_output:

  Logical. If `TRUE`, retain detailed ray, energy, and interception
  information, including the monthly ray discretization and interception
  matrices. If `FALSE`, only the main light-interception metrics are
  returned. Defaults to `FALSE`.

## Value

A list with two data frames, `sensors` and `cells`, with the same
columns as in the output of
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
(`light$sensors` and `light$cells`): `id_sensor` or `id_cell`, `e`,
`pacl` and `punobs`, and their direct and diffuse components if
`detailed_output = TRUE`. As in
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md),
energies are on a horizontal plane for sensors and on the slope plane
for cells. A data frame is empty if the geometry toward these targets
was not computed.

## Details

For each ray toward each target, the energy reaching the target is the
energy of the ray above the canopy times the transmittance of the
crossed crowns: \\\exp(-k \Omega \sum LAD \times length)\\ for turbid
media, or the product of the crown openness for porous envelopes. It is
null if a trunk intercepts the ray and `trunk_interception = TRUE`.

## See also

[precompute_sl_geometry](https://natheob.github.io/SamsaRaLight/reference/precompute_sl_geometry.md),
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)

## Examples

``` r
data_prenovel <- SamsaRaLight::data_prenovel

stand <- create_sl_stand(
  trees_inv = data_prenovel$trees,
  cell_size = 10,
  latitude = data_prenovel$info$latitude,
  slope = data_prenovel$info$slope,
  aspect = data_prenovel$info$aspect,
  north2x = data_prenovel$info$north2x
)
#> SamsaRaLight stand successfully created.

geometry <- precompute_sl_geometry(
  sl_stand = stand,
  monthly_radiations = data_prenovel$radiations
)
#> parallel mode: DISABLED | Threads = 1 (user) | Cores available = 4
#> SamsaRaLight geometry was computed successfully.

# Light on the cells for two crown LAD values
trees <- stand$trees

trees$crown_lad <- 0.5
head(estimate_pacl_from_geometry(geometry, trees)$cells)
#>   id_cell        e      pacl    punobs
#> 1       1 480.4950 0.1059478 0.1993924
#> 2       2 725.7215 0.1600196 0.4950759
#> 3       3 982.0738 0.2165446 0.7438652
#> 4       4 695.8185 0.1534261 0.5863522
#> 5       5 693.4853 0.1529116 0.2809002
#> 6       6 481.9791 0.1062751 0.5764287

trees$crown_lad <- 1
head(estimate_pacl_from_geometry(geometry, trees)$cells)
#>   id_cell        e       pacl    punobs
#> 1       1 269.1524 0.05934738 0.3559582
#> 2       2 516.1033 0.11379937 0.6961536
#> 3       3 830.7490 0.18317789 0.8793638
#> 4       4 516.0327 0.11378381 0.7906373
#> 5       5 451.4412 0.09954156 0.4315072
#> 6       6 348.3622 0.07681291 0.7975221
```
