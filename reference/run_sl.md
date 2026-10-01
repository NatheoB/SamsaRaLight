# Run SamsaRaLight radiative balance

Computes light interception and radiative balance for a forest stand
using the SamsaRaLight ray-tracing engine.

## Usage

``` r
run_sl(
  sl_stand,
  monthly_radiations,
  sensors_only = FALSE,
  detailed_output = FALSE,
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
  It contains tree, sensor, terrain, and grid geometry information.

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

- detailed_output:

  Logical. If `TRUE`, retain detailed ray, energy, and interception
  information, including the monthly ray discretization and interception
  matrices. If `FALSE`, only the main light-interception metrics are
  returned. Defaults to `FALSE`.

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

An object of class `"sl_output"` (a list) containing:

- `output`: A list containing the simulation results:

  - `light`: A list containing light-interception results for trees,
    ground cells, and sensors.

  - `monthly_rays`: The monthly ray discretization and associated
    radiation energies. Returned only when `detailed_output = TRUE`.

  - `interceptions`: Detailed interception matrices. Returned only when
    `detailed_output = TRUE`.

- `params`: A list containing the simulation parameters, including the
  simulated period, sky model, ray-discretization parameters, and
  interception-model parameters.

- `input`: A list containing the original `sl_stand` and
  `monthly_radiations` objects. Returned only when
  `include_input = TRUE`.

When `detailed_output = FALSE`, the main output tables contain:

- `light$sensors`: `id_sensor`, `e`, `pacl`, and `punobs`;

- `light$cells`: `id_cell`, `e`, `pacl`, and `punobs`;

- `light$trees`: `id_tree`, `epot`, `e`, `lci`, `eunobs`, and `rci`.

## Details

This function is the standard user interface to SamsaRaLight. Advanced
ray-tracing, sky-discretization, and interception parameters are
internally set to standard values and do not need to be provided.

Internally, `run_sl()` calls
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
using standard ray-tracing, sky-discretization, and interception
parameters:

- torus borders are enabled;

- crowns are represented as turbid media;

- the leaf extinction coefficient is `0.5`;

- the clumping factor is `1`;

- trunk interception is enabled;

- the minimum ray altitude angle is `10` degrees;

- direct rays use a starting offset of `0` degrees and an angular step
  of `5` degrees;

- diffuse rays use an angular step of `15` degrees;

- the Standard Overcast Sky model is used;

- the full year (days 1–365) is simulated.

These defaults are intended for standard SamsaRaLight simulations. Users
requiring control over these parameters can use
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
directly.

When `parallel_mode = TRUE`, SamsaRaLight uses OpenMP to parallelise
ray–target computations. BLAS is configured to use a single thread
during the simulation to avoid competing parallel execution.

## See also

[create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md),
[check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
[check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md),
[get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md),
[check_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/check_monthly_radiations.md),
[run_sl_advanced](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)

## Examples

``` r
data_prenovel <- SamsaRaLight::data_prenovel

stand <- create_sl_stand(
  trees = data_prenovel$trees,
  sensors = data_prenovel$sensors,
  cell_size = 10,
  latitude = data_prenovel$info$latitude,
  slope = data_prenovel$info$slope,
  aspect = data_prenovel$info$aspect,
  north2x = data_prenovel$info$north2x
)
#> SamsaRaLight stand successfully created.

out <- run_sl(
  sl_stand = stand,
  monthly_radiations = data_prenovel$radiations
)
#> parallel mode disabled because OpenMP was not available
#> SamsaRaLight simulation was run successfully.

str(out)
#> List of 3
#>  $ output:List of 1
#>   ..$ light:List of 3
#>   .. ..$ sensors:'data.frame':   0 obs. of  4 variables:
#>   .. .. ..$ id_sensor: int(0) 
#>   .. .. ..$ e        : num(0) 
#>   .. .. ..$ pacl     : num(0) 
#>   .. .. ..$ punobs   : num(0) 
#>   .. ..$ trees  :'data.frame':   333 obs. of  6 variables:
#>   .. .. ..$ id_tree: int [1:333] 242 258 272 116 157 37 89 224 241 92 ...
#>   .. .. ..$ epot   : num [1:333] 414830 684177 25792 411076 763878 ...
#>   .. .. ..$ e      : num [1:333] 176655 338207 1528 123157 359147 ...
#>   .. .. ..$ lci    : num [1:333] 0.574 0.506 0.941 0.7 0.53 ...
#>   .. .. ..$ eunobs : num [1:333] 150381 290852 615 97494 299592 ...
#>   .. .. ..$ rci    : num [1:333] 0.689 0.622 0.976 0.846 0.673 ...
#>   .. ..$ cells  :'data.frame':   100 obs. of  4 variables:
#>   .. .. ..$ id_cell: int [1:100] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. .. ..$ e      : num [1:100] 324 553 872 564 528 ...
#>   .. .. ..$ pacl   : num [1:100] 0.0714 0.122 0.1922 0.1243 0.1164 ...
#>   .. .. ..$ punobs : num [1:100] 0.238 0.616 0.848 0.723 0.369 ...
#>  $ input :List of 2
#>   ..$ sl_stand          :List of 7
#>   .. ..$ trees       :'data.frame':  333 obs. of  21 variables:
#>   .. .. ..$ id_tree       : int [1:333] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. .. ..$ species       : chr [1:333] "Abies alba" "Abies alba" "Abies alba" "Abies alba" ...
#>   .. .. ..$ x             : num [1:333] 77.8 62.7 84.1 59.1 33.4 ...
#>   .. .. ..$ y             : num [1:333] 71.5 65.6 95.6 98.3 39.8 ...
#>   .. .. ..$ z             : num [1:333] 7.51 6.89 10.05 10.33 4.18 ...
#>   .. .. ..$ dbh_cm        : num [1:333] 22.9 18.2 22.8 19.2 19.9 ...
#>   .. .. ..$ crown_type    : chr [1:333] "P" "P" "P" "P" ...
#>   .. .. ..$ h_m           : num [1:333] 14.8 13 15.8 11.7 14.7 ...
#>   .. .. ..$ hbase_m       : num [1:333] 3.31 4.74 4.91 4.21 3.83 ...
#>   .. .. ..$ hmax_m        : num [1:333] 3.31 4.74 4.91 4.21 3.83 ...
#>   .. .. ..$ rn_m          : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ re_m          : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rs_m          : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rw_m          : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ crown_lad     : num [1:333] 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 ...
#>   .. .. ..$ rxmax_m       : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rxmin_m       : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rymax_m       : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rymin_m       : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ added_to_fill : logi [1:333] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   .. .. ..$ crown_openness: num [1:333] NA NA NA NA NA NA NA NA NA NA ...
#>   .. ..$ sensors     : NULL
#>   .. ..$ cells       :'data.frame':  100 obs. of  6 variables:
#>   .. .. ..$ id_cell : num [1:100] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. .. ..$ x_center: num [1:100] 5 15 25 35 45 55 65 75 85 95 ...
#>   .. .. ..$ y_center: num [1:100] 95 95 95 95 95 95 95 95 95 95 ...
#>   .. .. ..$ z_center: num [1:100] 9.98 9.98 9.98 9.98 9.98 ...
#>   .. .. ..$ row     : int [1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ col     : int [1:100] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..$ core_polygon:List of 3
#>   .. .. ..$ df            :'data.frame': 15 obs. of  2 variables:
#>   .. .. .. ..$ x: num [1:15] 0.163 0.477 0.195 0.116 4.285 ...
#>   .. .. .. ..$ y: num [1:15] 7.62 14.6 84.6 91.04 96.43 ...
#>   .. .. ..$ sf            :Classes ‘sf’ and 'data.frame':    1 obs. of  2 variables:
#>   .. .. .. ..$ id      : int 1
#>   .. .. .. ..$ geometry:sfc_POLYGON of length 1; first list element: List of 1
#>   .. .. .. .. ..$ : num [1:16, 1:2] 0.163 0.477 0.195 0.116 4.285 ...
#>   .. .. .. .. ..- attr(*, "class")= chr [1:3] "XY" "POLYGON" "sfg"
#>   .. .. .. ..- attr(*, "sf_column")= chr "geometry"
#>   .. .. ..$ modify_polygon: chr "none"
#>   .. ..$ transform   :List of 10
#>   .. .. ..$ core_area_ha   : num 0.964
#>   .. .. ..$ core_batot_m2ha: num 36.7
#>   .. .. ..$ fill_around    : logi FALSE
#>   .. .. ..$ n_added_tree   : num 0
#>   .. .. ..$ new_area_ha    : num 1
#>   .. .. ..$ new_batot_m2ha : num 35.4
#>   .. .. ..$ epsg           : NULL
#>   .. .. ..$ shift_x        : num 0.0998
#>   .. .. ..$ shift_y        : num 0.38
#>   .. .. ..$ rotation_ccw   : num 0
#>   .. ..$ geometry    :List of 7
#>   .. .. ..$ cell_size: num 10
#>   .. .. ..$ n_cells_x: num 10
#>   .. .. ..$ n_cells_y: num 10
#>   .. .. ..$ latitude : num 46.5
#>   .. .. ..$ slope    : num 6
#>   .. .. ..$ aspect   : num 144
#>   .. .. ..$ north2x  : num 54
#>   .. ..$ inventory   :'data.frame':  333 obs. of  14 variables:
#>   .. .. ..$ id_tree   : int [1:333] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. .. ..$ species   : chr [1:333] "Abies alba" "Abies alba" "Abies alba" "Abies alba" ...
#>   .. .. ..$ x         : num [1:333] 77.7 62.6 84 59 33.3 ...
#>   .. .. ..$ y         : num [1:333] 71.1 65.2 95.2 97.9 39.4 ...
#>   .. .. ..$ dbh_cm    : num [1:333] 22.9 18.2 22.8 19.2 19.9 ...
#>   .. .. ..$ crown_type: chr [1:333] "P" "P" "P" "P" ...
#>   .. .. ..$ h_m       : num [1:333] 14.8 13 15.8 11.7 14.7 ...
#>   .. .. ..$ hbase_m   : num [1:333] 3.31 4.74 4.91 4.21 3.83 ...
#>   .. .. ..$ hmax_m    : logi [1:333] NA NA NA NA NA NA ...
#>   .. .. ..$ rn_m      : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ re_m      : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rs_m      : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ rw_m      : num [1:333] 3.02 3.09 2.83 2.52 2.82 ...
#>   .. .. ..$ crown_lad : num [1:333] 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 0.767 ...
#>   .. ..- attr(*, "class")= chr [1:2] "sl_stand" "list"
#>   ..$ monthly_radiations:'data.frame':   12 obs. of  3 variables:
#>   .. ..$ month  : int [1:12] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. ..$ Hrad   : num [1:12] 137 206 354 466 535 ...
#>   .. ..$ DGratio: num [1:12] 0.581 0.507 0.493 0.49 0.509 ...
#>  $ params:List of 13
#>   ..$ detailed_output   : logi FALSE
#>   ..$ start_day         : num 1
#>   ..$ end_day           : num 365
#>   ..$ soc               : logi TRUE
#>   ..$ use_torus         : logi TRUE
#>   ..$ turbid_medium     : logi TRUE
#>   ..$ extinction_coef   : num 0.5
#>   ..$ clumping_factor   : num 1
#>   ..$ trunk_interception: logi TRUE
#>   ..$ height_anglemin   : num 10
#>   ..$ direct_startoffset: num 0
#>   ..$ direct_anglestep  : num 5
#>   ..$ diffuse_anglestep : num 15
#>  - attr(*, "class")= chr [1:2] "sl_output" "list"
```
