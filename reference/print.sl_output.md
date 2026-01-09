# Print a `sl_output` object

Prints a concise overview of a SamsaRaLight simulation result, including
stand size, grid geometry, model configuration, and which outputs
(trees, cells, sensors, rays, interceptions) are available.

## Usage

``` r
# S3 method for class 'sl_output'
print(x, ...)
```

## Arguments

- x:

  An object of class `sl_output`.

- ...:

  Further arguments passed to or from other methods (currently ignored).

## Value

The object `x`, invisibly.

## Details

This method is designed to give a quick human-readable summary of what
was simulated and what was produced, without printing large data tables.

The printed output includes:

- Stand size (area, number of trees, cells, sensors)

- Grid geometry (number of cells and resolution)

- Radiation model parameters (latitude, time period, turbid medium,
  etc.)

- Which light and interception outputs are available

Use
[`summary.sl_output`](https://natheob.github.io/SamsaRaLight/reference/summary.sl_output.md)
to obtain numerical summaries of light, interception and energy
variables, and
[`plot.sl_output`](https://natheob.github.io/SamsaRaLight/reference/plot.sl_output.md)
for graphical inspection.

## See also

[`summary.sl_output`](https://natheob.github.io/SamsaRaLight/reference/summary.sl_output.md),
[`print.sl_stand`](https://natheob.github.io/SamsaRaLight/reference/print.sl_stand.md)
