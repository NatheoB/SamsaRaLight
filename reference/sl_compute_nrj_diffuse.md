# Compute energy of diffuse ray

Calculating SamsaraLight rays diffuse Energy in MJ/m2 of a plane
perpendicular to beam ray direction for a Standard Overcast Sky or
Uniform Overcast Sky are possible

## Usage

``` r
sl_compute_nrj_diffuse(
  soc,
  total_diffuse,
  height_angle_rad,
  diffuse_anglestep_rad
)
```

## Arguments

- soc:

  boolean - Standard Overcast Sky, if false: Uniform Overcast Sky

- total_diffuse:

  double - Total diffuse energy on a horizontal plan (in kWh.m-2)

- height_angle_rad:

  double - Angle between beam and soil (in radians)

- diffuse_anglestep_rad:

  double - Hour angle between two diffuse beams (in radians)

## Value

Energy per square meter of a horizontal plan (in MJ.m-2)
