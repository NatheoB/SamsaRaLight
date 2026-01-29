# Check coordinate columns and determine whether conversion is required

This function checks whether a data frame contains valid planar
coordinates (`x`, `y`) or geographic coordinates (`lon`, `lat`). If
geographic coordinates are detected, the user is informed that they must
be converted to a planar coordinate system using
[`create_xy_from_lonlat()`](https://natheob.github.io/SamsaRaLight/reference/create_xy_from_lonlat.md).

## Usage

``` r
check_coordinates(df, verbose = TRUE)
```

## Arguments

- df:

  A data.frame containing spatial coordinates.

- verbose:

  Logical. If TRUE, informative messages are printed.

## Value

Logical. Invisibly returns `TRUE` if coordinates must be converted from
longitude/latitude to planar coordinates, and `FALSE` if planar
coordinates are already present.

## Examples

``` r
df_xy <- data.frame(x = 1:3, y = 4:6)
check_coordinates(df_xy)
#> Planar coordinates successfully validated.

df_lonlat <- data.frame(lon = c(10, 11), lat = c(45, 46))
check_coordinates(df_lonlat, verbose = TRUE)
#> Geographic coordinates (`lon`, `lat`) detected.
```
