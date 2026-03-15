# Plot a SamsaRaLight output

Visualize ground light and tree energy metrics for a `sl_output` object.

## Usage

``` r
# S3 method for class 'sl_output'
plot(
  x,
  ...,
  what_trees = c("compet", "intercepted", "potential"),
  what_cells = c("relative", "absolute"),
  show_trees = TRUE,
  direct_energy = NULL
)
```

## Arguments

- x:

  An object of class `sl_output`, returned by
  [`run_sl()`](https://natheob.github.io/SamsaRaLight/reference/run_sl.md).

- ...:

  Additional arguments passed to lower-level plotting functions.

- what_trees:

  Character; which tree metric to plot. Choices are:

  "compet"

  :   Light competition index (LCI), reversed viridis scale.

  "intercepted"

  :   Intercepted energy (MJ).

  "potential"

  :   Potential intercepted energy (MJ).

  Default is "compet".

- what_cells:

  Character; which cell (ground) metric to plot. Choices are:

  "relative"

  :   Proportion of above canopy light (PACL)).

  "absolute"

  :   Energy on the ground (MJ/m2).

  Default is "relative".

- show_trees:

  Logical; whether to display trees on top of the ground light map.
  Default is TRUE.

- direct_energy:

  Logical or NULL. If NULL (default), total radiation outputs are
  plotted (direct + diffuse). If TRUE, only direct radiation components
  are plotted. If FALSE, only diffuse radiation components are plotted.
  This option requires `detailed_output = TRUE` when running the
  simulation

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(sl_output_object)
plot(sl_output_object, what_trees = "potential", what_cells = "absolute", show_trees = FALSE)
} # }
```
