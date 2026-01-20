# Plot a from-above view of a tree inventory

Visualizes a forest inventory as ellipses representing tree crowns.

## Usage

``` r
plot_inventory(trees_inv, show_id = TRUE)
```

## Arguments

- trees_inv:

  A data.frame of trees that passed
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md).

- show_id:

  Logical; if TRUE (default), displays tree identifiers at crown
  centers.

## Value

A ggplot object displaying the trees as ellipses in a from-above view.
