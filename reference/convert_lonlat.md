# Convert data.frame from lon/lat to UTM planar coordinates

Convert data.frame from lon/lat to UTM planar coordinates

## Usage

``` r
convert_lonlat(df, epsg_utm)
```

## Arguments

- df:

  A data.frame containing `lon` and `lat` columns

- epsg_utm:

  EPSG code of the UTM projection

## Value

Data.frame with added `x` and `y` columns (meters)
