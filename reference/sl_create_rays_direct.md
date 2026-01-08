# Create direct rays

Create SamsaraLight diffuse rays of a plane with a classical sky
hemisphere divided by meridians and parallels can use Standard Overcast
Sky or Uniform Overcast Sky

## Usage

``` r
sl_create_rays_direct(
  latitude_rad,
  declination_rad,
  nrj_direct_months,
  height_anglemin_rad,
  direct_anglestep_rad,
  slope_rad,
  bottom_azimut_rad,
  southazimut_ccw_rad,
  direct_startoffset_rad
)
```

## Arguments

- latitude_rad:

  double - Latitude of the plot (in radians)

- declination_rad:

  double - Declination angle in radians:

- nrj_direct_months:

  12 double vect - Monthly direct energy on a horizontal plan (in
  kWh.m-2)

- height_anglemin_rad:

  double - Angle minimum between beam and soil (in radians)

- direct_anglestep_rad:

  double - Hour angle between two diffuse beams (in radians)

- slope_rad:

  double - Slope of the plan (in radians)

- bottom_azimut_rad:

  double - Azimut of the vector orthogonal to the ground

- southazimut_ccw_rad:

  double - Azimuth of south counterclockwise from x axis in the (x,y)
  system (in radians)

- direct_startoffset_rad:

  double - Angle at which to start first direct ray (in radians)

## Value

list of 3 elements : horizontal energy (double), slope energy (double)
and diffuse rays (data.frame)
