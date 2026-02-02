# Rotate 2D vectors by an angle (counter-clockwise)

Internal utility to rotate 2D coordinates by a given angle (in radians),
using a standard rotation matrix.

## Usage

``` r
rotate_vec_ccw(x, y, theta)
```

## Arguments

- x:

  Numeric vector of x coordinates.

- y:

  Numeric vector of y coordinates.

- theta:

  Rotation angle in radians (counter-clockwise).

## Value

A data.frame with rotated coordinates:

- x:

  Rotated x coordinates

- y:

  Rotated y coordinates

## Details

The rotation follows: \$\$ x' = x \cos(\theta) - y \sin(\theta) \$\$
\$\$ y' = x \sin(\theta) + y \cos(\theta) \$\$
