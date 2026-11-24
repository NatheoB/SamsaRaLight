# Precompute the geometry of light interception in a stand

Runs, once for a given stand, the part of the SamsaRaLight ray tracing
that does not depend on the crown leaf area density (LAD) or openness:
for each ray toward each target (sensor or cell), the crowns crossed by
the ray, the path length across each crown, and whether a trunk
intercepts the ray. The light reaching the targets can then be estimated
very quickly for any crown LAD or openness with
[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md),
e.g. when calibrating crown LAD from light measurements, where the same
stand is run many times with different LAD values.

## Usage

``` r
precompute_sl_geometry(
  sl_stand,
  monthly_radiations,
  sensors_only = FALSE,
  use_torus = TRUE,
  height_anglemin = 10,
  direct_startoffset = 0,
  direct_anglestep = 5,
  diffuse_anglestep = 15,
  soc = TRUE,
  start_day = 1,
  end_day = 365,
  parallel_mode = FALSE,
  n_threads = NULL,
  verbose = TRUE
)
```

## Arguments

- sl_stand:

  An object of class `"sl_stand"` describing the forest stand, typically
  created with
  [create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md).
  It contains tree, sensor, terrain, and grid geometry information (see
  [validate_sl_stand](https://natheob.github.io/SamsaRaLight/reference/validate_sl_stand.md)).

- monthly_radiations:

  A data frame containing monthly horizontal radiation (`Hrad`, in MJ
  m\\^{-2}\\) and diffuse-to-global radiation ratio (`DGratio`),
  typically obtained with
  [get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md).
  The input is validated internally using
  [check_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/check_monthly_radiations.md).

- sensors_only:

  Logical. If `TRUE`, compute light interception only for sensors. This
  can substantially reduce computation time when only sensor-level
  results are required. Defaults to `FALSE`.

- use_torus:

  Logical. If `TRUE`, use a torus system to handle stand borders.
  Defaults to `TRUE`.

- height_anglemin:

  Numeric. Minimum altitude angle of rays, in degrees. Defaults to `10`.

- direct_startoffset:

  Numeric. Starting angle of the first direct ray, in degrees. Defaults
  to `0`.

- direct_anglestep:

  Numeric. Angular step between direct rays, in degrees. Defaults to
  `5`.

- diffuse_anglestep:

  Numeric. Angular step between diffuse rays, in degrees. Defaults to
  `15`.

- soc:

  Logical. If `TRUE`, use the Standard Overcast Sky model; if `FALSE`,
  use the Uniform Overcast Sky model. Defaults to `TRUE`.

- start_day:

  Numeric. First day of the simulated vegetative period, between 1
  and 365. Defaults to `1`.

- end_day:

  Numeric. Last day of the simulated vegetative period, between 1
  and 365. Must be greater than or equal to `start_day`. Defaults to
  `365`.

- parallel_mode:

  Logical. If `TRUE`, ray–target computations are parallelised using
  OpenMP. If `FALSE`, the model runs in single-thread mode. BLAS is
  automatically forced to single-thread mode during the simulation to
  avoid competition between BLAS and OpenMP. Parallel execution can
  substantially reduce computation time for large stands or fine ray
  discretization. Defaults to `FALSE`.

- n_threads:

  Integer or `NULL`. Number of CPU threads to use when
  `parallel_mode = TRUE`. If `NULL`, OpenMP automatically selects the
  number of available threads. If supplied, must be a positive integer.
  Defaults to `NULL`.

- verbose:

  Logical. If `TRUE`, print informative messages during the simulation,
  including OpenMP status. Defaults to `TRUE`.

## Value

An object of class `"sl_geometry"` (a list), to be used in
[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md),
and that can be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). It contains:

- `sensors`, `cells`: The geometry toward the sensors and the cells
  (`NULL` if there are no sensors, or if `sensors_only = TRUE` for the
  cells). Each one is a list with:

  - `targets`: One row per target (`id_sensor` or `id_cell`), with the
    energy of the rays reaching the target without crossing any crown or
    trunk (`e_unobs_direct`, `e_unobs_diffuse`, in MJ/m2 on a horizontal
    plane for sensors and on the slope plane for cells, as in
    [run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)).

  - `pairs`: One row per ray x target pair intercepted at least once by
    a crown or a trunk, with the index of the target (`target`, row in
    `targets`) and of the ray (`ray`, row in `rays`), whether a trunk
    intercepts the ray (`blocked`), and the number of crowns crossed
    (`n_trees`).

  - `crossings`: One row per crown crossed by the ray of a pair, with
    the index of the pair (`pair`, row in `pairs`), of the tree (`tree`,
    row in `sl_stand$trees`) and the path length across the crown
    (`length`, in m, summed over crown parts).

- `rays`: One row per ray, with its energy above the canopy (in MJ/m2,
  on a horizontal plane `e_horizontal` and on the slope plane `e_slope`)
  and its type (`direct`).

- `id_tree`: Identifiers of the trees, in the order of `sl_stand$trees`
  (order of the crown LAD or openness to give in
  [estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md)).

- `energies`: Energies above the canopy, as computed for
  [run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md).

- `params`: The settings used to compute the geometry.

The object should not be modified: filter the output of
[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md)
instead.

## Details

The energy reaching a target does not depend on the order in which a ray
crosses the crowns: it is the energy of the ray above the canopy times
the product of the transmittances of the crossed crowns (Beer-Lambert
law for a turbid medium, crown openness for a porous envelope), and is
null if a trunk intercepts the ray. The light on sensors and cells
estimated by
[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md)
is thus the same as with
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
(up to rounding errors) for the same stand, radiations and ray settings.
The light intercepted by trees (`e`, `epot`, `lci`, ...) depends on the
order of interceptions and cannot be estimated from the geometry.

The geometry depends on the stand (tree positions, dimensions and crown
shapes, sensors, cells), the radiations and the ray settings, but not on
the crown LAD and openness, the extinction coefficient, the clumping
factor or the interception by trunks, which are chosen in
[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md).
The geometry must be computed again if the stand, the radiations or the
ray settings change. Trunks are always included in the geometry.

The size of the geometry increases with the number of targets, rays and
crowns crossed by each ray: use `sensors_only = TRUE` if the light on
the cells is not needed.

## See also

[estimate_pacl_from_geometry](https://natheob.github.io/SamsaRaLight/reference/estimate_pacl_from_geometry.md),
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

str(geometry, max.level = 2)
#> List of 6
#>  $ sensors : NULL
#>  $ cells   :List of 3
#>   ..$ targets  :'data.frame':    100 obs. of  3 variables:
#>   ..$ pairs    :'data.frame':    46314 obs. of  4 variables:
#>   ..$ crossings:'data.frame':    204871 obs. of  3 variables:
#>  $ rays    :'data.frame':    491 obs. of  3 variables:
#>   ..$ e_horizontal: num [1:491] 0.959 1.343 1.699 2.026 2.321 ...
#>   ..$ e_slope     : num [1:491] 1.82 2.21 2.57 2.89 3.18 ...
#>   ..$ direct      : logi [1:491] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ id_tree : int [1:333] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ energies: Named num [1:4] 2440 2095 2289 2061
#>   ..- attr(*, "names")= chr [1:4] "slope_direct" "slope_diffuse" "horizontal_direct" "horizontal_diffuse"
#>  $ params  :List of 9
#>   ..$ sensors_only      : logi FALSE
#>   ..$ start_day         : num 1
#>   ..$ end_day           : num 365
#>   ..$ soc               : logi TRUE
#>   ..$ use_torus         : logi TRUE
#>   ..$ height_anglemin   : num 10
#>   ..$ direct_startoffset: num 0
#>   ..$ direct_anglestep  : num 5
#>   ..$ diffuse_anglestep : num 15
#>  - attr(*, "class")= chr [1:2] "sl_geometry" "list"
```
