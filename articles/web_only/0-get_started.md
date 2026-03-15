# Get started

*SamsaRaLight* is documented through a series of progressive tutorials
that guide users from simple simulations to advanced applications.
Although it is advisable to follow the tutorials in the correct order,
they are still independent of each other.

1.  [A first minimal
    case](https://natheob.github.io/SamsaRaLight/articles/web_only/1-minimal_case.md) -
    *Compute tree light interception by symmetric crowns in an
    axis-aligned rectangle plot*

2.  [Understand ray
    discretization](https://natheob.github.io/SamsaRaLight/articles/web_only/2-discretisation_directdiffuse.md) -
    *Discretization of direct and diffuse radiations into rays*

3.  *Represent tree crowns with asymmetric shapes and handle more
    complex inventory geometries.*

    1.  [Non-axis-aligned rectangular stand from GPS
        data](https://natheob.github.io/SamsaRaLight/articles/web_only/3a-complex_cases_IRRES.md):
        Convert geographic coordinates and create an axis-aligned
        virtual stand.

    2.  [Asymmetric crowns in a non-axis-aligned
        stand](https://natheob.github.io/SamsaRaLight/articles/web_only/3b-complex_cases_bechefa.md):
        Use asymmetric crown geometries and fill the surrounding area
        with virtual trees.

    3.  [Non-rectangular inventory
        zone](https://natheob.github.io/SamsaRaLight/articles/web_only/3c-complex_cases_cloture.md):
        Work with complex inventory polygons.

4.  [Understand ray
    transmission](https://natheob.github.io/SamsaRaLight/articles/web_only/4-transmission_interception.md) -
    *Explore transmission models for light propagation through tree
    crowns, including porous envelope and turbid medium approaches.*

5.  *Estimate the light reaching virtual sensors placed in the stand.*

    1.  [Implement virtual
        sensors](https://natheob.github.io/SamsaRaLight/articles/web_only/5a-virtual_sensors_concept.md):
        Introduce the concept and implementation of virtual sensors in
        the plot.

    2.  [Calibrate the LAD
        parameter](https://natheob.github.io/SamsaRaLight/articles/web_only/5b-virtual_sensors_lad.md):
        Use virtual sensor measurements to calibrate the leaf area
        density parameter

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
