# Validate an sl_geometry object

Checks the structure of an object created by
[precompute_sl_geometry](https://natheob.github.io/SamsaRaLight/reference/precompute_sl_geometry.md)
(format version 2): its components, the links between `targets`,
`pairs`, `rays` and `trees`, and the dimensions and content of the
matrices `L`, `C` and `W`, for the pairs not intercepted and intercepted
by a trunk.

## Usage

``` r
validate_sl_geometry(x)
```

## Arguments

- x:

  An object of class `"sl_geometry"`.

## Value

Invisibly `TRUE`, or an error.
