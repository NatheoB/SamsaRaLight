# Example forest inventory datasets for SamsaRaLight

These datasets provide example forest inventories used for light
interception simulations with the SamsaRaLight package. Each dataset is
a named list with 5 elements: `trees`, `sensors`, `core_polygon`,
`radiations`, and `info`.

## Usage

``` r
data(data_prenovel)
data(data_IRRES1)
data(data_bechefa)
data(data_cloture20)

data_prenovel

data_IRRES1

data_bechefa

data_cloture20
```

## Format

A named list with 5 elements:

- trees:

  A `data.frame` with tree-level data:

  id_tree

  :   Unique id of the tree (integer).

  species

  :   Latin species name (character).

  x, y

  :   Coordinates of the base of the tree in meters (numeric).

  dbh_cm

  :   Diameter at breast height in cm (numeric).

  crown_type

  :   Crown shape type, e.g. "P", "E", or complex types like "8E"
      (character).

  h_m

  :   Height of the tree trunk in meters (numeric).

  hbase_m

  :   Height of crown base in meters (numeric).

  hmax_m

  :   Height at which crown radius is maximum, NA if not used (numeric).

  rn_m, rs_m, re_m, rw_m

  :   Crown radii in meters (numeric).

  crown_lad

  :   Leaf area density of the crown (m2/m3) (numeric).

- sensors:

  A `data.frame` with sensor data:

  id_sensor

  :   Unique sensor ID (integer).

  x, y

  :   Coordinates of the sensor in meters (numeric).

  h_m

  :   Height of the sensor in meters (numeric).

  pacl, pacl_direct, pacl_diffuse

  :   Proportion of above-canopy light measured at the sensor (numeric).

  May be `NULL` if no sensors are present.

- core_polygon:

  A `data.frame` with vertices of the inventory polygon:

  x, y

  :   Coordinates of polygon vertices in meters (numeric).

- radiations:

  A `data.frame` of monthly radiation:

  month

  :   Month number (integer).

  Hrad

  :   Monthly radiation in MJ (numeric).

  DGratio

  :   Diffuse-to-global radiation ratio (numeric).

- info:

  A named list with site information:

  latitude, longitude

  :   Coordinates of the site in decimal degrees (numeric).

  slope

  :   Mean slope in degrees (numeric).

  aspect

  :   Aspect in degrees from north (numeric).

  north2x

  :   Angle from north to x-axis clockwise in degrees (numeric).

An object of class `list` of length 5.

An object of class `list` of length 5.

An object of class `list` of length 5.

An object of class `list` of length 5.

## Source

Courbaud Benoit (INRAe LESSEM Grenoble)

Gauthier Ligot (Gembloux Agro-Bio Tech)

Gauthier Ligot (Gembloux Agro-Bio Tech)

Gauthier Ligot (Gembloux Agro-Bio Tech)
