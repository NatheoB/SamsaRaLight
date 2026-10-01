# Compute advanced SamsaRaLight radiative balance

Runs the full light interception and radiative balance simulation for a
virtual forest stand with advanced ray-tracing and interception
parameters. This function provides full control over ray discretization,
sky model, crown interception model, trunk interception, and parallel
computation.

## Usage

``` r
run_sl_advanced(
  sl_stand,
  monthly_radiations,
  sensors_only = FALSE,
  use_torus = TRUE,
  turbid_medium = TRUE,
  extinction_coef = 0.5,
  clumping_factor = 1,
  trunk_interception = TRUE,
  height_anglemin = 10,
  direct_startoffset = 0,
  direct_anglestep = 5,
  diffuse_anglestep = 15,
  soc = TRUE,
  start_day = 1,
  end_day = 365,
  detailed_output = FALSE,
  parallel_mode = FALSE,
  n_threads = NULL,
  verbose = TRUE
)
```

## Arguments

- sl_stand:

  An object of class `"sl_stand"` representing the virtual forest stand.
  Each row of `sl_stand$trees` represents a tree, with required and
  optional columns describing crown geometry, tree height, crown radius,
  crown openness, leaf area density, and related attributes. See
  [validate_sl_stand](https://natheob.github.io/SamsaRaLight/reference/validate_sl_stand.md).

- monthly_radiations:

  A data frame containing monthly horizontal radiation (`Hrad`) and
  diffuse-to-global radiation ratios (`DGratio`), typically computed
  with
  [get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md).

- sensors_only:

  Logical. If `TRUE`, compute light interception only for sensors.
  Defaults to `FALSE`.

- use_torus:

  Logical. If `TRUE`, use a torus system to handle stand borders.
  Defaults to `TRUE`.

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

- detailed_output:

  Logical. If `TRUE`, retain detailed ray, energy, and interception
  information in the output. If `FALSE`, only the main
  light-interception metrics are retained. Defaults to `FALSE`.

- parallel_mode:

  Logical. If `TRUE`, ray–target computations are parallelised using
  OpenMP. If `FALSE`, the model runs in single-thread mode. Defaults to
  `FALSE`.

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

  - `light`: Light-interception results for trees, cells, and sensors.
    When `detailed_output = FALSE`, these contain the main output
    metrics only.

  - `monthly_rays`: The generated monthly ray discretization and
    associated radiation energies. Returned only when
    `detailed_output = TRUE`.

  - `interceptions`: Detailed tree/cell interception matrices. Returned
    only when `detailed_output = TRUE`.

- `params`: A list containing the simulation parameters, including the
  simulated period, sky model, ray discretization, crown interception
  model, and interception parameters.

- `input`: A list containing the original `sl_stand` and
  `monthly_radiations` objects. Returned only when
  `include_input = TRUE`.

When `detailed_output = FALSE`, the main output tables contain:

- `light$sensors`: `id_sensor`, `e`, `pacl`, and `punobs`;

- `light$cells`: `id_cell`, `e`, `pacl`, and `punobs`;

- `light$trees`: `id_tree`, `epot`, `e`, `lci`, `eunobs`, and `rci`.

## Details

For typical use, see
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md),
which provides standard values for the ray-discretization parameters.

This advanced function exposes all ray-tracing parameters used by
[`create_sl_rays`](https://natheob.github.io/SamsaRaLight/reference/create_sl_rays.md)
and the interception model used by the underlying C++ simulation. It is
intended for users who need fine control over ray discretization, sky
conditions, crown representation, or computational settings.

Before running the simulation, the function validates the stand, monthly
radiation data, interception-model configuration, logical and numeric
parameters, simulation period, and number of OpenMP threads.

BLAS and OpenMP thread counts are configured to avoid competing parallel
execution. When `parallel_mode = TRUE`, the number of OpenMP threads can
be controlled with `n_threads`.

For most users,
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md) is
recommended because it provides standard ray-discretization parameters.
