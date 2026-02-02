# Plot a from-above view of a tree inventory

Visualizes a forest inventory as ellipses representing tree crowns.

## Usage

``` r
plot_inventory(trees_inv, core_polygon_df = NULL, show_id = TRUE)
```

## Arguments

- trees_inv:

  A data.frame of trees that passed
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md).

- core_polygon_df:

  Optional data.frame defining the core inventory polygon. Must contain
  columns `x` and `y`.

- show_id:

  Logical; if TRUE (default), displays tree identifiers at crown
  centers.

## Value

A ggplot object displaying the trees in a from-above view.

## Details

Because the north2x variable is unknown, trees are plotted as circles by
considering the mean radius on the four cardinals.
