# Validate a SamsaRaLight simulation output object

Checks the internal consistency of an object returned by
`run_sl_advanced`.

## Usage

``` r
validate_sl_output(x)
```

## Arguments

- x:

  Object of class `"sl_output"`.

## Value

Invisibly returns TRUE if validation passes, stops with informative
message otherwise.

## Details

Validates that:

- Object inherits from class `"sl_output"`.

- `light` component contains `trees`, `cells`, and `sensors` as
  data.frames.

- All essential columns (`id_tree`, `id_cell`, `id_sensor`, energy
  variables) exist.
