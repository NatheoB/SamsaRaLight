# Run SamsaRaLight radiative balance

This function computes light interception and radiative balance for a
forest stand using the SamsaRaLight ray-tracing engine.

## Usage

``` r
run_sl(
  sl_stand,
  monthly_radiations,
  latitude,
  sensors_only = FALSE,
  use_torus = TRUE,
  turbid_medium = TRUE,
  detailed_output = FALSE,
  parallel_mode = FALSE,
  n_threads = NULL,
  verbose = TRUE
)
```

## Arguments

- sl_stand:

  An object of class `"sl_stand"` describing the forest stand, created
  with
  [create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md).
  It contains trees, sensors, terrain, and grid geometry.

- monthly_radiations:

  A data.frame of monthly horizontal radiation (`Hrad`, in MJ
  m\\^{-2}\\) and diffuse-to-global ratio (`DGratio`), typically
  obtained using
  [get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md)
  and checked using
  [check_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/check_monthly_radiations.md)..

- latitude:

  Numeric. Latitude of the stand (degrees).

- sensors_only:

  Logical. If `TRUE`, compute light interception only for sensors (much
  faster).

- use_torus:

  Logical. If `TRUE`, stand borders are treated using a torus (periodic)
  geometry, mimicking an infinite forest, representative of the virtual
  stand. If `FALSE`, borders are open (surrounded by grassland).

- turbid_medium:

  Logical. If `TRUE`, tree crowns are treated as a turbid medium
  (Beer–Lambert law using `crown_lad`). If `FALSE`, crowns are treated
  as porous envelopes (using `crown_openness`).

- detailed_output:

  Logical. If `TRUE`, the output contains detailed diffuse/direct
  energies in the `light` datasets, full interception matrices
  `interceptions` and output of ray discretization `monthy_rays`. If
  `FALSE`, only total energies are returned (recommended for most uses).

- parallel_mode:

  logical. If TRUE, ray–target computations are parallelised using
  OpenMP. If FALSE, the model runs in single-thread mode. SamsaRaLight
  uses OpenMP for ray–target parallelisation. To avoid competition
  between OpenMP and BLAS (matrix algebra libraries), BLAS is
  automatically forced to single-thread mode during the simulation.
  Using `parallel_mode = TRUE` is strongly recommended for large stands
  or fine ray discretisation, as computation time scales almost linearly
  with the number of available CPU cores.

- n_threads:

  integer or NULL. Number of CPU threads to use when
  `parallel_mode = TRUE`. If NULL (default), OpenMP automatically
  selects the number of available cores. If provided, must be a positive
  integer.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

An object of class `"sl_output"`, containing:

- light:

  A list of data.frames with simulated light interception:

  - `trees`: light intercepted by trees

  - `cells`: light received by ground cells

  - `sensors`: light received by sensors

- info:

  A list of metadata about the simulation (latitude, sky type, torus
  use, etc.).

- monthly_rays:

  (only if `detailed_output = TRUE`) Discretization of monthly
  radiations

- interceptions:

  (only if `detailed_output = TRUE`) interception matrices between trees
  and rays for each cell/sensor

## Details

It is the **standard user interface** of SamsaRaLight. Advanced
ray-tracing and sky discretization parameters are internally set to
robust defaults and do not need to be provided.

Internally, `run_sl()` calls the advanced engine
[`run_sl_advanced()`](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
with fixed ray-tracing and sky discretization.

You should normally **not** use `SamsaRaLight:::run_sl_advanced()`
directly unless you are developing new ray-tracing configurations or
doing methodological work.

## See also

[create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md),
[check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
[check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md),
[get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md),
[check_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/check_monthly_radiations.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data_prenovel <- SamsaRaLight::data_prenovel

stand <- create_sl_stand(
  trees = data_prenovel$trees,
  sensors = data_prenovel$sensors,
  cell_size = 5,
  slope = 10,
  aspect = 180,
  north2x = 0
)

rad <- get_monthly_radiations(45.8, 3.1)

out <- run_sl(
  sl_stand = stand,
  monthly_radiations = rad,
  latitude = 45.8
)

out$light$trees
} # }
```
