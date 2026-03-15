# Check the format and validity of a sensor position data.frame

This function checks whether a sensor data.frame is correctly formatted
to be used as input for the ray-tracing model. It verifies the presence,
type, and validity of mandatory variables describing the position and
height of sensors within a forest stand.

## Usage

``` r
check_sensors(sensors, verbose = TRUE)
```

## Arguments

- sensors:

  A data.frame with one row per sensor and the following columns:

  id_sensor

  :   Unique identifier of the sensor (numeric or character, no
      duplicates)

  x

  :   X position of the sensor (numeric, meters)

  y

  :   Y position of the sensor (numeric, meters)

  h_m

  :   Height above ground of the sensor (numeric, meters)

- verbose:

  Logical; if TRUE, informative messages are printed.

## Value

Invisibly returns TRUE if all checks pass.

## Examples

``` r
if (FALSE) { # \dontrun{
sensors <- data.frame(
  id_sensor = 1:3,
  x = c(10, 20, 30),
  y = c(5, 15, 25),
  h_m = c(1.5, 2.0, 1.8)
)
check_sensors(sensors)
} # }
```
