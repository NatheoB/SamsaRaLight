# Create planar (x, y) coordinates from longitude / latitude

This function converts geographic coordinates (`lon`, `lat`) into planar
coordinates (`x`, `y`) using an automatically selected UTM projection.

## Usage

``` r
create_xy_from_lonlat(df)
```

## Arguments

- df:

  A data.frame containing geographic coordinates (`lon` and `lat`).

## Value

A list with the following elements:

- `df`:

  The input data.frame with added `x` and `y` columns (meters, UTM).

- `epsg`:

  EPSG code of the UTM projection used.

## Details

The input data.frame must contain geographic coordinates (`lon` and
`lat`). Planar coordinates (`x`, `y`) are automatically computed using
the UTM zone inferred from the mean longitude and hemisphere inferred
from the mean latitude.
