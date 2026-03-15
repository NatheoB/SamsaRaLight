# Validate a SamsaRaLight stand object

This function checks the internal consistency and structure of an object
of class `"sl_stand"`, as returned by
[create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md).
It verifies that all required components are present and correctly
formatted, that the embedded tree inventory conforms to the rules
enforced by
[check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
and that all trees and sensors are within the stand limits.

## Usage

``` r
validate_sl_stand(x)
```

## Arguments

- x:

  An object expected to be of class `"sl_stand"`.

## Value

Invisibly returns `TRUE` if the stand is valid.

## Details

The following validations are performed:

- The object inherits from class `"sl_stand"`.

- The top-level components `trees`, `sensors`, `cells`, `core_polygon`,
  `transform`, `geometry` and `inventory` are present.

- The `trees` data.frame passes
  [check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md).

- The `sensors` data.frame passes
  [check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md).

- The `cells` data.frame contains columns `x_center`, `y_center`,
  `z_center`, and `id_cell`.

- The `geometry` list contains `cell_size`, `n_cells_x`, `n_cells_y`,
  `slope`, `aspect`, and `north2x`.

- All trees and sensors lie within the bounds of the rectangular stand.

## See also

[create_sl_stand](https://natheob.github.io/SamsaRaLight/reference/create_sl_stand.md),
[check_inventory](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md),
[check_sensors](https://natheob.github.io/SamsaRaLight/reference/check_sensors.md),
[run_sl](https://natheob.github.io/SamsaRaLight/reference/run_sl.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data_prenovel <- SamsaRaLight::data_prenovel
stand <- create_sl_stand(data_prenovel$trees, cell_size = 5,
                         slope = 10, aspect = 180, north2x = 0)

validate_sl_stand(stand)  # returns TRUE invisibly
} # }
```
