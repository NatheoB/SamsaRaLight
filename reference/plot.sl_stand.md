# Plot a SamsaRaLight virtual stand

This function plots a virtual forest stand (`sl_stand`) produced by
SamsaRaLight. It can display a top-down view with tree crowns or a
side/top view with cells and trees.

## Usage

``` r
# S3 method for class 'sl_stand'
plot(x, top_down = FALSE, transparency = TRUE)
```

## Arguments

- x:

  An object of class `sl_stand`.

- top_down:

  Logical, if TRUE, creates a top-down view with multiple directions
  (south, north, west, east).

- transparency:

  Logical, if TRUE, trees are semi-transparent.

## Value

A `ggplot` object representing the stand.
