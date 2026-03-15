# Validate a SamsaRaLight simulation output object

Performs structural and internal consistency checks on an object
returned by
[`run_sl_advanced()`](https://natheob.github.io/SamsaRaLight/reference/run_sl_advanced.md)
or
[`run_sl()`](https://natheob.github.io/SamsaRaLight/reference/run_sl.md).

## Usage

``` r
validate_sl_output(x)
```

## Arguments

- x:

  Object expected to inherit from class `"sl_output"`.

## Value

Invisibly returns `TRUE` if validation passes. Stops with an informative
error message otherwise.

## Details

This function is called internally at the end of a simulation to ensure
that the returned object is valid and complete before being provided to
the user.

The function checks that:

- The object inherits from class `"sl_output"`.

- Top-level components `$output` and `$input` exist.

- `$output$light` contains `trees`, `cells`, and `sensors`.

- These elements are data.frames.

- Required identifier and energy columns are present.

- Energy-related columns are numeric.

If detailed outputs were requested, the presence and structure of
`$output$monthly_rays` and `$output$interceptions` are also verified
when available.
