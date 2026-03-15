# Validate monthly radiation input

Checks that a monthly radiation data.frame is correctly formatted and
physically valid for use by the light interception and ray-tracing
model. The table must contain exactly 12 months of radiation data.

## Usage

``` r
check_monthly_radiations(x, verbose = TRUE)
```

## Arguments

- x:

  A data.frame with monthly radiation values, typically produced by
  [`get_monthly_radiations`](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md).

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

Invisibly returns `TRUE` if all checks pass.

## Details

The input must contain the following columns:

- month:

  Integer month number (1-12)

- Hrad:

  Monthly global horizontal irradiation (MJ/m2)

- DGratio:

  Diffuse-to-global radiation ratio (unitless, 0-1)

The function checks:

- 1:

  The object is a data.frame.

- 2:

  Required columns are present.

- 3:

  There are exactly 12 months.

- 4:

  Each month (1-12) is present exactly once.

- 5:

  Data are numeric and finite.

- 6:

  Hrad strictly positive.

- 7:

  DGratio between 0 and 1.

- 8:

  Months are in increasing order.
