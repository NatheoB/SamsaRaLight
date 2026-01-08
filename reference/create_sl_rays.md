# Compute direct and diffuse rays of a growing season

Compute direct and diffuse rays of a growing season

## Usage

``` r
create_sl_rays(
  monthly_rad,
  latitude,
  start_day = 1,
  end_day = 365,
  soc = TRUE,
  slope = 0,
  north_to_x_cw = 90,
  aspect = 0,
  height_anglemin = 10,
  direct_startoffset = 0,
  direct_anglestep = 5,
  diffuse_anglestep = 15
)
```

## Arguments

- monthly_rad:

  data.frame - Monthly horizontal radiation (Hrad) and diffuse to global
  ratio (DGratio) Computed with function samsaRa::sl_get_monthlyrad()

- latitude:

  double - Latitude of the plot (in degrees)

- start_day:

  integer between 1 and 365 - First day of the vegetative period

- end_day:

  integer between 1 and 365 - Last day of the vegetative period

- soc:

  boolean - Standard Overcast Sky, if false: Uniform Overcast Sky

- slope:

  double - Slope of the plot (in degrees)

- north_to_x_cw:

  double - Angle from North to x axis clockwise. (in degrees) Default
  correspond to a Y axis oriented toward the North.

- aspect:

  double - Angle of slope bottom on the compass from the North,
  clockwise rotation (in degrees) northern aspect : 0, eastern aspect :
  90, southern aspect : 180, western aspect : 270

- height_anglemin:

  double - Angle minimum between beam and soil (in degrees)

- direct_startoffset:

  double - Angle at which to start first direct ray (in degrees)

- direct_anglestep:

  double - Hour angle between two direct beams (in degrees)

- diffuse_anglestep:

  double - Hour angle between two diffuse beams (in degrees)

## Value

list of 3 elements : horizontal energy (double), slope energy (double)
and rays (data.frame) with n rows and 5 columns:

- azimutAzimut of the ray in radians

- height_angleAngle between beam and soil (in radians)

- eEnergy of ray before crossing the canopy (in MJ.m-2)

- directtrue if the ray is direct false if it is diffuse
