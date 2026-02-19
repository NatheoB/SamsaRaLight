# Get started

*SamsaRaLight* is documented through a series of progressive tutorials
that guide users from simple simulations to advanced applications.
Although it is advisable to follow the tutorials in the correct order,
they are still independent of each other.

1.  [A first minimal
    case](https://natheob.github.io/SamsaRaLight/articles/1-minimal_case.md) -
    *Compute tree light interception by symmetric crowns in an
    axis-aligned rectangle plot*
2.  [Understand ray
    discretization](https://natheob.github.io/SamsaRaLight/articles/2-discretisation_directdiffuse.md) -
    *Discretization of direct and diffuse radiations into rays*
3.  [More complex crown shapes and
    inventories](https://natheob.github.io/SamsaRaLight/articles/3-complex_cases.md) -
    *Represent tree crowns with asymmetric shapes and use
    non-axis-aligned GPS inventory or non-rectangular inventory zone*
4.  [Understand ray
    transmission](https://natheob.github.io/SamsaRaLight/articles/4-transmission_interception.md) -
    *Transmission models : porous envelope and turbid medium crowns*
5.  [Virtual
    sensors](https://natheob.github.io/SamsaRaLight/articles/5-virtual_sensors.md) -
    *Estimate light arriving towards a virtual sensor*

Below is a typical example workflow.

``` r
library(SamsaRaLight)

# Load example inventory
data("data_cloture20")

# Create virtual stand
stand <- create_sl_stand(
  trees_inv = data_cloture20$trees,
  cell_size = 5,
  latitude = 50.03617,
  slope = 0,
  aspect = 0,
  north2x = 90,
  sensors = data_cloture20$sensors,
  core_polygon_df = data_cloture20$core_polygon,
  fill_around = TRUE
)

# Observe input stand
summary(stand)
plot(stand)
plot(stand, top_down = TRUE)

# Get radiation data from PVGIS online database
rad <- get_monthly_radiations(
  latitude = 50.03617, 
  longitude = 5.20634
)

# Run the model
out <- run_sl(
  sl_stand = stand,
  monthly_rad = rad
)

# Plot results
summary(out)
plot(out)
plot(out, show_trees = FALSE)
```
