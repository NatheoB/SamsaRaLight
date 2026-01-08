# Information for IRRES1 inventory

Information about IRRES1 inventory in Ardennes (Belgium) used as an
example for informing the inventory zone with a core polygon.

## Usage

``` r
data(data_IRRES1)
```

## Format

`data_IRRES1`, a named list with 5 elements:

- trees:

  data.frame that contains for each tree, the information about its
  location, species, size, size of the crown and crown information.

  id_tree

  :   Unique id of the tree. (integer)

  species

  :   Species Latin name. (character)

  x, y

  :   Coordinates of the base of the tree in the forest stand in m.
      (double)

  dbh_cm

  :   Diameter at breast height (1.30m) of the trunk of the tree in cm.
      (double)

  crown_type

  :   Type of the crown between paraboloid (P) and ellispoid (E).
      (character)

  h_m

  :   Height of the tree trunk in m. (double)

  hbase_m

  :   Crown base height of the tree (i.e. height at wich the crown
      start) in m. (double)

  hmax_m

  :   Height at which the crown radius is maximum. Set it to NA in case
      of simple crown shapes ("P" and "E"), otherwise, your values won't
      be considered as hmax is automatically computed.

  rn_m, rs_m, re_m, rw_m

  :   Biggest radius of the tree crown in m (double)

  crown_openness

  :   Crown Openness of the tree (no unit) (i.e. Fraction of the energy
      of a light ray crossing the crown that is intercepted). Used when
      computing interception with a crown considered as a porous
      envelop. (double)

  crown_lad

  :   Leaf Area Density of the tree crown in m2/m3 (i.e. surface of
      leave per volume of crown, considering an homogeneous crown). Used
      when computing interception with a crown considered as a turbid
      medium. (double)

- species_color:

  named vector that contains color HEX code for each species

- sensors:

  data.frame that contains information about the stand sensors

- core_polygon:

  data.frame containing vertices of the tree inventory zone, described
  by the coordinates (x, y) of each vertex edge.

- radiation:

  data.frame that contains monthly radiations.

  month

  :   Unique id of the tree. (integer)

  Hrad

  :   Monthly radiation in a horizontal plane in MJ. (double)

  DGratio

  :   Ratio between monthly diffuse and global energies. (double)

- info:

  A named numeric vector with site information:

  latitude

  :   Latitude of the site (decimal degrees).

  longitude

  :   Longitude of the site (decimal degrees).

  size_x

  :   Horizontal extent of the site (meters). Here, it is set to NA as
      the inventory is defined by the core polygon.

  size_y

  :   Vertical extent of the site (meters). Here, it is set to NA as the
      inventory is defined by the core polygon.

  slope

  :   Mean slope of the terrain (degrees).

  aspect

  :   Aspect (degrees from north).

  north_to_x_cw

  :   Angle from north to the x-axis (clockwise, degrees).

## Source

Gauthier Ligot (Gembloux Agro-Bio Tech)
