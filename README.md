
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
n_cells <- 20 # Number of cells within one length
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
# trees <- as.data.table(SamsaRaLight::data_trees_prenovel %>% 
#                          dplyr::mutate(rn_m = r_m,
#                                        re_m = r_m,
#                                        rs_m = r_m,
#                                        rw_m = r_m,
#                                        crown_type = "8E",
#                                        hmax_m = hbase_m + 1/2*(h_m - hbase_m)) %>% 
#                          dplyr::select(-r_m)) # Try to give a data.table to the function to avoid converting inside the function trees
trees <- as.data.table(SamsaRaLight::data_trees_bechefa)
trees
#>      id_tree     species     x    y    dbh_cm crown_type  h_m hbase_m hmax_m
#>   1:     103       abies 163.1 67.3  68.43663         8E 38.5    16.6   28.2
#>   2:     615 pseudotsuga  66.8 41.4  95.49297         8E 50.2    14.0   33.3
#>   3:     102 pseudotsuga 159.2 58.2 111.72677         8E 48.2    10.8   27.1
#>   4:     708 pseudotsuga  43.8 34.2 116.81973         8E 51.0    15.8   27.0
#>   5:     707 pseudotsuga  51.4 34.1  99.94930         8E 50.5    16.0   26.7
#>  ---                                                                        
#> 197:      15       fagus 195.0 79.1  16.55211         8E 13.7     0.6    2.3
#> 198:     703       fagus  42.2 48.7  28.64789         8E 21.1     1.0    2.1
#> 199:     704       fagus  55.8 44.5  20.05352         8E 11.2     1.8    1.8
#> 200:     511       fagus  75.3 57.7  21.00845         8E  8.4     1.3    1.8
#> 201:     700       fagus  20.2 76.5  44.88169         8E 18.4     1.0    1.6
#>      rn_m rs_m re_m rw_m crown_openess crown_lad
#>   1: 4.32 4.12 3.70 5.20         0.035       0.6
#>   2: 7.01 6.45 5.87 4.15         0.035       0.6
#>   3: 8.38 6.72 4.69 6.51         0.035       0.6
#>   4: 7.37 7.43 9.85 5.44         0.035       0.6
#>   5: 6.73 5.16 3.40 5.76         0.035       0.6
#>  ---                                            
#> 197: 3.27 4.80 3.00 4.65         0.035       0.6
#> 198: 6.27 4.40 6.85 6.02         0.035       0.6
#> 199: 3.51 5.65 4.85 0.64         0.035       0.6
#> 200: 5.12 2.60 3.60 3.72         0.035       0.6
#> 201: 2.03 4.07 2.22 6.04         0.035       0.6
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

``` r
microbenchmark::microbenchmark(
  "sl" = SamsaRaLight::sl_run(
    trees, monthly_rad,
    latitude = latitude, slope = slope, 
    aspect = aspect, north_to_x_cw = north_to_x_cw,
    cell_size = cell_size, n_cells = n_cells,
    use_rcpp = T,
    turbid_medium = TRUE
  ))
#> Unit: milliseconds
#>  expr     min       lq     mean   median      uq     max neval
#>    sl 44.6723 47.10285 50.38295 48.85395 53.0911 60.5104   100
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
#>     id_tree           epot              e           
#>  Min.   :  1.0   Min.   :     0   Min.   :     0.0  
#>  1st Qu.:210.0   1st Qu.:  3678   1st Qu.:   308.6  
#>  Median :413.0   Median :  9284   Median :  1934.2  
#>  Mean   :403.3   Mean   : 26908   Mean   :  8973.6  
#>  3rd Qu.:607.0   3rd Qu.: 26233   3rd Qu.:  6891.9  
#>  Max.   :815.0   Max.   :274064   Max.   :152579.9
```

`cells`: Light coming to each cell of the plot

- `e`: Total energy arriving to the cell (in MJ/year)

- `erel`: Relative energy coming to the cell compared to the one above
  canopy. It ranges from 0 (no light on the ground) to 1 (no energy has
  been intercepted by above trees).

``` r
summary(out$cells)
#>     id_cell            e               erel        
#>  Min.   :  1.0   Min.   :-246.2   Min.   :-0.5362  
#>  1st Qu.:100.8   1st Qu.: 459.2   1st Qu.: 1.0000  
#>  Median :200.5   Median : 459.2   Median : 1.0000  
#>  Mean   :200.5   Mean   : 414.3   Mean   : 0.9023  
#>  3rd Qu.:300.2   3rd Qu.: 459.2   3rd Qu.: 1.0000  
#>  Max.   :400.0   Max.   : 459.2   Max.   : 1.0000
```
