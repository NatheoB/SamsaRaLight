
# SamsaRaLight

The goal of SamsaRaLight is to estimate the light intercepted by a tree
and the light on the ground in a given forest stand.

It uses the ray tracing model SamsaraLight (Courbaud et al. 2003, 2015)
to compute the direct and diffuse rays within a given period of a year,
at a given latitude. Then, it simulates the ray coming from the sky
towards cells of the plot and it computes the interception of all the
rays by each consecutive tree crowns with given dimensions.

The SamsaraLight model was initially implemented within the Java
platform Capsis (<https://capsis.cirad.fr/capsis/help_en/samsaralight>).
However, for the sake of fast computing and easier use within R
paradigm, we needed to reformulate the model as an R package. We used
both `data.table` R package and Rcpp with C++ script for fast computing.

## Installation

You can install the development version of SamsaRaLight from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("NatheoB/SamsaRaLight")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SamsaRaLight)
library(data.table)
```

### Set the stand geometry

You first need to precise the geometry of the stand with:

``` r
# Coordinates of the stand (for monthly radiation and rays geometry)
latitude <- 46 # The latitude of the stand (Y-coord in WGS85)
longitude <- 2 # The longitud eof the plot (X-coord in WGS85)

# The slope of the stand
slope <- 6

# Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
# northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
aspect <- 144

# Angle from North to x axis clockwise. (in degrees)
# Default correspond to a Y axis oriented toward the North.
north_to_x_cw <- 54

# Considering a squared plot, number and size of the cells composing the grid
n_cells <-  10 # Number of cells within one length
cell_size <- 10 # Size of the length of a cell
```

### Create the trees dataset

To be able to run SamsaRaLight, you need to provide a well formatted
data.frame containing precise information about each tree of the plot.

The trees data.frame should contain all the below variables, with the
correct name and type of the column:

- `id_tree`: Unique id of the tree. (integer)

- `species`: Species Latin name. (character)

- `x`, `y`, `z`: Coordinates of the base of the tree in the forest stand
  in m. (double)

- `dbh_cm`: Diameter at breast height (1.30m) of the trunk of the tree
  in cm. (double)

- `height_m`: Height of the tree trunk in m. (double)

- `cbh_m`: Crown base height of the tree (i.e. height at wich the crown
  start) in m. (double)

- `cradius_m`: Biggest radius of the tree crown in m (double)

- `crown_type`: Type of the crown between “paraboloid” (1) and
  “ellispoid (2). (integer)

- `crown_lad`: Leaf Area Density of the tree crown in m2/m3 (i.e.
  surface of leave per volume of crown, considering an homogeneous
  crown). Used when computing interception with a crown considered as a
  turbid medium. (double)

- `crown_p`: Crown Openness of the tree (no unit) (i.e. Fraction of the
  energy of a light ray crossing the crown that is intercepted). Used
  when computing interception with a crown considered as a porous
  envelop. (double)

You can find an example tree dataset
`SamsaRaLight::data_trees_prenovel`. It is an uneven-aged mountain stand
of fir, spruce, beech located in Prenovel (Jura, France).

``` r
trees <- as.data.table(SamsaRaLight::data_trees_prenovel) # Try to give a data.table to the function to avoid converting inside the function
trees
#>      id_tree     species       x       y       z  dbh_cm height_m  cbh_m
#>   1:       1  Abies alba 77.7336 71.0808  7.4709 22.9320  14.8120 3.3073
#>   2:       2  Abies alba 62.5783 65.1864  6.8514 18.2397  12.9612 4.7429
#>   3:       3  Abies alba 84.0060 95.2483 10.0110 22.8472  15.8187 4.9055
#>   4:       4  Abies alba 58.9521 97.9011 10.2898 19.1539  11.7033 4.2050
#>   5:       5  Abies alba 33.3423 39.4053  4.1416 19.8852  14.6597 3.8346
#>  ---                                                                    
#> 329:     329  Abies alba 83.6769 29.9210  3.1448 16.1720  12.9273 2.7850
#> 330:     330  Abies alba 20.1081 22.6300  2.3785 13.9037  10.4369 4.6049
#> 331:     331  Abies alba 90.7703 27.8787  2.9302  9.4128   8.1643 2.3558
#> 332:     332 Picea abies 90.1975 15.2572  1.6036 17.1754  16.0397 6.7220
#> 333:     333  Abies alba 94.8614 14.0905  1.4810 15.3389   9.3093 2.3800
#>      cradius_m crown_type crown_lad crown_p
#>   1:    3.0204          1       0.5     0.2
#>   2:    3.0896          1       0.5     0.2
#>   3:    2.8350          1       0.5     0.2
#>   4:    2.5196          1       0.5     0.2
#>   5:    2.8208          1       0.5     0.2
#>  ---                                       
#> 329:    2.7448          1       0.5     0.2
#> 330:    2.3478          1       0.5     0.2
#> 331:    1.8538          1       0.5     0.2
#> 332:    2.3005          1       0.5     0.2
#> 333:    2.6004          1       0.5     0.2
```

### Get monthly radiation

You will also need to provide a dataset with two radiation values for
each of the 12 months:

- `month`: the month between 1 and 12

- `Hrad`: monthly global radiation on a horizontal plane in MJ

- `DGratio`: ratio of monthly diffuse energy to global energy (needed to
  compute direct energy)

You can fetch this dataset for a given latitude and longitude using the
function `get_monthly_rad()`. It fetches monthly data between given
start and end year (range between 2005 and 2020) from the PVGIS database
<https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system_en>.
Here, you can choose to average the monthly values between all the years
between 2005 and 2020 or to choose values of a given year.

``` r
monthly_rad <- SamsaRaLight::get_monthly_rad(latitude, longitude)
monthly_rad
#> # A tibble: 12 × 3
#>    month  Hrad DGratio
#>    <int> <dbl>   <dbl>
#>  1     1  134.   0.619
#>  2     2  198.   0.563
#>  3     3  356.   0.506
#>  4     4  475.   0.495
#>  5     5  556.   0.497
#>  6     6  621.   0.476
#>  7     7  669.   0.433
#>  8     8  587.   0.428
#>  9     9  439.   0.441
#> 10    10  281.   0.509
#> 11    11  161.   0.568
#> 12    12  123.   0.591
```

### Run SamsaRaLight

Now, given the stand geometry, the trees dataset and the monthly
radiation, you can run SamsaRaLight using the function `sl_run()`.

``` r
out <- SamsaRaLight::sl_run(
  trees, monthly_rad,
  latitude = latitude, slope = slope, 
  aspect = aspect, north_to_x_cw = north_to_x_cw,
  cell_size = cell_size, n_cells = n_cells
)
```

The function returns a list with two dataframes:

`trees`: Light interception for each tree

- `epot`: Potential energy intercepted by the tree without its
  neighbours (in MJ/year)

- `e`: Energy intercepted by the tree when considering competition with
  neighbours (in MJ/year)

- `lci`: Light competition index, computed as $lci = 1-\frac{e}{epot}$.
  It ranges from 0 to 1 with 0 being a tree without competition for
  light (all the energy that the tree can intercept from coming rays are
  intercepted) and 1 being a tree totally in competition for light (all
  the energy from the rays coming to the tree crown has been intercepted
  by the neighbours).

``` r
summary(out$trees)
#>     id_tree         epot              e                 lci         
#>  Min.   :  1   Min.   :  8080   Min.   :   681.3   Min.   :0.06347  
#>  1st Qu.: 84   1st Qu.:138419   1st Qu.: 24617.7   1st Qu.:0.53748  
#>  Median :167   Median :250658   Median : 71870.9   Median :0.69330  
#>  Mean   :167   Mean   :280214   Mean   :114401.8   Mean   :0.66048  
#>  3rd Qu.:250   3rd Qu.:397578   3rd Qu.:179545.3   3rd Qu.:0.81312  
#>  Max.   :333   Max.   :934841   Max.   :695745.0   Max.   :0.97177
```

`cells`: Light coming to each cell of the plot

- `e`: Total energy arriving to the cell (in MJ/year)

- `erel`: Relative energy coming to the cell compared to the one above
  canopy. It ranges from 0 (no light on the ground) to 1 (no energy has
  been intercepted by above trees).

``` r
 summary(out$cells)
#>     id_cell             e               erel        
#>  Min.   :  1.00   Min.   : 253.8   Min.   :0.05497  
#>  1st Qu.: 25.75   1st Qu.: 558.3   1st Qu.:0.12092  
#>  Median : 50.50   Median : 749.6   Median :0.16235  
#>  Mean   : 50.50   Mean   : 828.3   Mean   :0.17941  
#>  3rd Qu.: 75.25   3rd Qu.:1008.0   3rd Qu.:0.21832  
#>  Max.   :100.00   Max.   :1695.5   Max.   :0.36723
```

## Other parameters

Cf turbid medium, use_rcpp, use_torus…

## Speed comparison

### With 10m cells

``` r
microbenchmark::microbenchmark(
  "r" = SamsaRaLight::sl_run(
    trees, monthly_rad,
    latitude = latitude, slope = slope, 
    aspect = aspect, north_to_x_cw = north_to_x_cw,
    cell_size = 10, n_cells = 10,
    use_rcpp = FALSE
  ),
  "rcpp" = SamsaRaLight::sl_run(
    trees, monthly_rad,
    latitude = latitude, slope = slope, 
    aspect = aspect, north_to_x_cw = north_to_x_cw,
    cell_size = 10, n_cells = 10,
    use_rcpp = TRUE
  ),
  times = 10
)
#> Unit: milliseconds
#>  expr        min         lq       mean     median         uq        max neval
#>     r 10301.5682 10909.2967 11698.6835 11298.7528 12003.7309 15066.7850    10
#>  rcpp   165.7329   174.3573   187.8694   185.2281   206.9962   212.2717    10
#>  cld
#>   a 
#>    b
```

### With 5m cells

``` r
# microbenchmark::microbenchmark(
#   "r" = SamsaRaLight::sl_run(
#     trees, monthly_rad,
#     latitude = latitude, slope = slope, 
#     aspect = aspect, north_to_x_cw = north_to_x_cw,
#     cell_size = 5, n_cells = 20,
#     use_rcpp = FALSE
#   ),
#   "rcpp" = SamsaRaLight::sl_run(
#     trees, monthly_rad,
#     latitude = latitude, slope = slope, 
#     aspect = aspect, north_to_x_cw = north_to_x_cw,
#     cell_size = 5, n_cells = 20,
#     use_rcpp = TRUE
#   ),
#   times = 10
# )
```
