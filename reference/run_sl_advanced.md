# Compute advanced SamsaRaLight radiative balance

This function runs the full light interception and radiative balance
simulation for a virtual forest stand with advanced parameters. It
allows customization of ray discretization, sky type and trunk
interception.

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

  An object of class `"sl_stand"` representing the virtual stand. Each
  row is a tree with required and optional columns describing crown
  geometry, height, crown radius, crown openness, LAD, etc. See
  [validate_sl_stand](https://natheob.github.io/SamsaRaLight/reference/validate_sl_stand.md).

- monthly_radiations:

  data.frame of monthly horizontal radiation (Hrad) and diffuse to
  global ratio (DGratio), computed with
  [get_monthly_radiations](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md).

- sensors_only:

  logical, if TRUE, compute interception only for sensors

- use_torus:

  logical, if TRUE, use torus system for borders, else open grassland

- turbid_medium:

  logical, if TRUE, crowns are considered turbid medium, else porous
  envelope

- extinction_coef:

  Numeric scalar. Leaf extinction coefficient controlling the
  probability that a ray is intercepted by foliage. It represents the
  effective light attenuation per unit leaf area and is linked to
  average leaf orientation. Higher values increase interception (default
  = 0.5).

- clumping_factor:

  Numeric scalar controlling the aggregation of leaves within the crown
  volume. A value of 1 corresponds to a homogeneous (random) foliage
  distribution; values \< 1 indicate clumped foliage, and values \> 1
  indicate more regular spacing. This modifies effective light
  interception in the turbid medium model (default = 1).

- trunk_interception:

  logical, if TRUE, account for trunk interception

- height_anglemin:

  numeric, minimum altitude angle for rays (degrees)

- direct_startoffset:

  numeric, starting angle of first direct ray (degrees)

- direct_anglestep:

  numeric, hour angle step between direct rays (degrees)

- diffuse_anglestep:

  numeric, hour angle step between diffuse rays (degrees)

- soc:

  logical, if TRUE, use Standard Overcast Sky; if FALSE, Uniform
  Overcast Sky

- start_day:

  integer, first day of the vegetative period (1–365)

- end_day:

  integer, last day of the vegetative period (1–365)

- detailed_output:

  logical, if TRUE, include detailed rays, energies, and interception
  matrices

- parallel_mode:

  logical. If TRUE, ray–target computations are parallelised using
  OpenMP. If FALSE, the model runs in single-thread mode.

- n_threads:

  integer or NULL. Number of CPU threads to use when
  `parallel_mode = TRUE`. If NULL (default), OpenMP automatically
  selects the number of available cores. If provided, must be a positive
  integer.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

An object of class `"sl_output"` (list) containing:

- `light`: list with simulation outputs for trees, cells, and sensors

- `info`: list with run metadata (latitude, days, sky type, etc.)

- `monthly_rays` (if detailed_output = TRUE): ray discretization per
  month

- `interceptions` (if detailed_output = TRUE): tree/cell interception
  matrices

## Details

For typical use, see the simpler
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md)
wrapper that sets standard discretization parameters for most users.

This advanced function exposes all ray tracing parameters and is
intended for users who need full control over ray discretization and
modeling options. For most users, see
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md)
which wraps this function with default parameters suitable for standard
runs.
