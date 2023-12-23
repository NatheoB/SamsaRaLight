
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

## Small example

This is a basic example which shows you how to solve a common problem:

``` r
library(SamsaRaLight)
library(data.table)
library(dplyr)
```

### Set the stand geometry

You first need to precise the geometry of the stand with:

``` r
# Coordinates of the stand (for monthly radiation and rays geometry)
latitude <- 46 # The latitude of the stand (Y-coord in WGS85)
longitude <- 2 # The longitude of the plot (X-coord in WGS85)

# The slope of the stand
slope <- 6

# Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
# northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
aspect <- 144

# Angle from North to x axis clockwise. (in degrees)
# Default correspond to a Y axis oriented toward the North.
north_to_x_cw <- 54

# Considering a squared plot, number and size of the cells composing the grid
n_cells <- 10 # Number of cells within one length
cell_size <- 10 # Size of the length of a cell
```

### Create the trees dataset

To be able to run SamsaRaLight, you need to provide a well formatted
data.frame containing precise information about each tree of the plot.

The trees data.frame should contain all the below variables, with the
correct name and type of the column:

- `id_tree`: Unique id of the tree. (integer)

- `species`: Species Latin name. (character)

- `x`, `y`: Coordinates of the base of the tree in the forest stand
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
trees <- as.data.table(SamsaRaLight::data_trees_prenovel %>%
                         dplyr::mutate(rn_m = r_m,
                                       re_m = r_m,
                                       rs_m = r_m,
                                       rw_m = r_m,
                                       crown_type = "E",
                                       hmax_m = hbase_m + 1/2*(h_m - hbase_m)) %>%
                         dplyr::select(-r_m)) 
# Try to give a data.table to the function to avoid converting inside the function trees
# trees <- as.data.table(SamsaRaLight::data_trees_bechefa)
trees
#>      id_tree     species       x       y  dbh_cm crown_type     h_m hbase_m
#>   1:       1  Abies alba 77.7336 71.0808 22.9320          E 14.8120  3.3073
#>   2:       2  Abies alba 62.5783 65.1864 18.2397          E 12.9612  4.7429
#>   3:       3  Abies alba 84.0060 95.2483 22.8472          E 15.8187  4.9055
#>   4:       4  Abies alba 58.9521 97.9011 19.1539          E 11.7033  4.2050
#>   5:       5  Abies alba 33.3423 39.4053 19.8852          E 14.6597  3.8346
#>  ---                                                                       
#> 329:     329  Abies alba 83.6769 29.9210 16.1720          E 12.9273  2.7850
#> 330:     330  Abies alba 20.1081 22.6300 13.9037          E 10.4369  4.6049
#> 331:     331  Abies alba 90.7703 27.8787  9.4128          E  8.1643  2.3558
#> 332:     332 Picea abies 90.1975 15.2572 17.1754          E 16.0397  6.7220
#> 333:     333  Abies alba 94.8614 14.0905 15.3389          E  9.3093  2.3800
#>      crown_openess crown_lad   rn_m   re_m   rs_m   rw_m   hmax_m
#>   1:           0.2       0.5 3.0204 3.0204 3.0204 3.0204  9.05965
#>   2:           0.2       0.5 3.0896 3.0896 3.0896 3.0896  8.85205
#>   3:           0.2       0.5 2.8350 2.8350 2.8350 2.8350 10.36210
#>   4:           0.2       0.5 2.5196 2.5196 2.5196 2.5196  7.95415
#>   5:           0.2       0.5 2.8208 2.8208 2.8208 2.8208  9.24715
#>  ---                                                             
#> 329:           0.2       0.5 2.7448 2.7448 2.7448 2.7448  7.85615
#> 330:           0.2       0.5 2.3478 2.3478 2.3478 2.3478  7.52090
#> 331:           0.2       0.5 1.8538 1.8538 1.8538 1.8538  5.26005
#> 332:           0.2       0.5 2.3005 2.3005 2.3005 2.3005 11.38085
#> 333:           0.2       0.5 2.6004 2.6004 2.6004 2.6004  5.84465
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

Monthly radiation for Prenovel example is stored within the package
using `SamsaRaLight::data_rad_prenovel`. The dataset is fetched using
the above function.

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
    cell_size = cell_size, n_cells = n_cells,
    use_rcpp = T,
    turbid_medium = TRUE
  )
```

    #> Unit: milliseconds
    #>  expr      min      lq     mean  median       uq      max neval
    #>    sl 198.7417 207.083 234.1921 213.138 221.1911 464.6904   100

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
#>     id_tree         epot               e         
#>  Min.   :  1   Min.   :  10072   Min.   :   562  
#>  1st Qu.: 84   1st Qu.: 163738   1st Qu.: 22453  
#>  Median :167   Median : 299307   Median : 57445  
#>  Mean   :167   Mean   : 330761   Mean   :119931  
#>  3rd Qu.:250   3rd Qu.: 470347   3rd Qu.:171867  
#>  Max.   :333   Max.   :1085693   Max.   :794042
```

`cells`: Light coming to each cell of the plot

- `e`: Total energy arriving to the cell (in MJ/year)

- `erel`: Relative energy coming to the cell compared to the one above
  canopy. It ranges from 0 (no light on the ground) to 1 (no energy has
  been intercepted by above trees).

``` r
summary(out$cells)
#>     id_cell          x_center     y_center     z_center            e         
#>  Min.   :  1.00   Min.   : 5   Min.   : 5   Min.   :0.5255   Min.   : 195.7  
#>  1st Qu.: 25.75   1st Qu.:25   1st Qu.:25   1st Qu.:2.6276   1st Qu.: 451.9  
#>  Median : 50.50   Median :50   Median :50   Median :5.2552   Median : 575.0  
#>  Mean   : 50.50   Mean   :50   Mean   :50   Mean   :5.2552   Mean   : 645.2  
#>  3rd Qu.: 75.25   3rd Qu.:75   3rd Qu.:75   3rd Qu.:7.8828   3rd Qu.: 772.2  
#>  Max.   :100.00   Max.   :95   Max.   :95   Max.   :9.9849   Max.   :1329.8  
#>       erel        
#>  Min.   :0.04238  
#>  1st Qu.:0.09789  
#>  Median :0.12454  
#>  Mean   :0.13974  
#>  3rd Qu.:0.16724  
#>  Max.   :0.28803
```
