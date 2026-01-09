# Plot a from-above view of a tree inventory

Visualizes a forest inventory as ellipses representing tree crowns.
Taller trees are drawn above shorter ones. Optionally colors trees by
species.

## Usage

``` r
plot_inventory(trees_inv, transparency = TRUE, show_id = TRUE)
```

## Arguments

- trees_inv:

  A data.frame of trees that passed
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md).
  Must contain at least:

  `x`, `y`

  :   Tree coordinates (meters).

  `rn_m`, `rs_m`, `re_m`, `rw_m`

  :   Crown radii (meters).

  `h_m`

  :   Total height (meters) for plotting order.

  `id_tree`

  :   Optional if show_id is FALSE; Tree identifier used for labeling.

  `species`

  :   Optional; species name (character) for coloring.

- transparency:

  Logical; if TRUE (default), uses alpha = 0.6, otherwise alpha = 1
  (opaque).

- show_id:

  Logical; if TRUE (default), displays tree identifiers at crown
  centers.

## Value

A ggplot object displaying the trees as ellipses in a from-above view.

## Details

Time for plotting can be long because trees are plotted one by one in
order to respect height layering for both crowns and labels.
