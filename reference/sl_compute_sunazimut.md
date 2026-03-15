# Compute sun azimut

Computaion of sun azimut for a given height angle reference system with
angle origin on X \> 0 axis and trigonometric rotation

## Usage

``` r
sl_compute_sunazimut(
  latitude_rad,
  declination_rad,
  hour_angle_rad,
  height_angle_rad,
  southazimut_ccw_rad
)
```

## Arguments

- latitude_rad:

  double - Latitude of the plot (in radians)

- declination_rad:

  double - Declination angle in radians: angle between the equator and a
  line drawn from the centre of the Earth to the centre of the sun

- hour_angle_rad:

  double - Angle which the sun moves across the sky (in radians)

- height_angle_rad:

  double - Angle between beam and soil (in radians)

- southazimut_ccw_rad:

  double - Azimuth of south counterclockwise from x axis in the (x,y)
  system (in radians)
