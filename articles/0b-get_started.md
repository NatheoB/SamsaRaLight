# Get started

*SamsaRaLight* is documented through a series of progressive tutorials
that guide users from simple simulations to advanced applications.
Although it is advisable to follow the tutorials in the correct order,
they are still independent of each other.

1.  [A first minimal
    case](https://natheob.github.io/SamsaRaLight/articles/vignette(%221-minimal_case%22)) -
    *Compute tree light interception by symmetric crowns in an
    axis-aligned rectangle plot*
2.  [Understand ray
    discretization](https://natheob.github.io/SamsaRaLight/articles/vignette(%222-discretisation_directdiffuse%22)) -
    *Discretization of direct and diffuse radiations into rays*
3.  [More complex crown shapes and
    inventories](https://natheob.github.io/SamsaRaLight/articles/vignette(%223-complex_cases%22)) -
    *Represent tree crowns with asymmetric shapes and use
    non-axis-aligned GPS inventory or non-rectangular inventory zone*
4.  [Understand ray
    transmission](https://natheob.github.io/SamsaRaLight/articles/vignette(%224-transmission_interception%22)) -
    *Transmission models : porous envelope and turbid medium crowns*
5.  [Virtual
    sensors](https://natheob.github.io/SamsaRaLight/articles/vignette(%225-virtual_sensors%22)) -
    *Estimate light arriving towards a virtual sensor*

Below is a minimal example workflow.

``` r
library(SamsaRaLight)

# Load example inventory
data("data_prenovel")

# Create virtual stand
stand <- create_sl_stand(
  trees_inv = data_prenovel$trees,
  cell_size = 5,
  latitude = 46.52666,
  slope = 6,
  aspect = 144,
  north2x = 54,
)

# Observe input stand
summary(stand)
plot(stand)

# Get radiation data from PVGIS online database
rad <- get_monthly_radiations(
  latitude = 46.52666, 
  longitude = 5.82765
)

# Run the model
out <- run_sl(
  sl_stand = stand,
  monthly_rad = rad
)

# Plot results
summary(out)
plot(out)
```
