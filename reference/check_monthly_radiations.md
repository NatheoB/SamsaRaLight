# Validate monthly radiation input

Checks that a monthly radiation data.frame is correctly formatted and
physically valid for use by the light interception and ray-tracing
model. The table must contain exactly 12 months of radiation data.

## Usage

``` r
check_monthly_radiations(x)
```

## Arguments

- x:

  A data.frame with monthly radiation values, typically produced by
  [`get_monthly_radiations`](https://natheob.github.io/SamsaRaLight/reference/get_monthly_radiations.md).

## Value

Invisibly returns `TRUE` if all checks pass.

## Details

The input must contain the following columns:

- monthInteger month number (1–12)

- HradMonthly global horizontal irradiation (MJ m⁻²)

- DGratioDiffuse-to-global radiation ratio (unitless, 0–1)

The function checks:

- The object is a data.frame

- Required columns are present

- There are exactly 12 months

- Each month (1–12) is present exactly once

- Data are numeric and finite

- Hrad ≥ 0

- 0 ≤ DGratio ≤ 1

- Months are in increasing order
