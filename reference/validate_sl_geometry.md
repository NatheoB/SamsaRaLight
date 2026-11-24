# Validate an sl_geometry object

Checks the structure of an object created by
[precompute_sl_geometry](https://natheob.github.io/SamsaRaLight/reference/precompute_sl_geometry.md),
and the links between its tables (indices of targets, rays, pairs and
trees, and number of crowns crossed per pair).

## Usage

``` r
validate_sl_geometry(x)
```

## Arguments

- x:

  An object of class `"sl_geometry"`.

## Value

Invisibly `TRUE`, or an error.
