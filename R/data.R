#' Information for Prenovel inventory
#'
#' @description Information about Prenovel marteloscope in Jura (France) used 
#' as an example for light interception by symetrical crowns.
#' 
#' @format `data_prenovel`, a named list with 3 elements:
#'  \describe{
#'    \item{trees}{data.frame that contains for each tree, the information about 
#'      its location, species, size, size of the crown and crown information.
#'      \describe{
#'        \item{id_tree}{Unique id of the tree. (integer)}
#'        \item{species}{Species Latin name. (character)}
#'        \item{x, y}{Coordinates of the base of the tree in the forest stand in m. (double)}
#'        \item{dbh_cm}{Diameter at breast height (1.30m) of the trunk of the tree in cm. (double)}
#'        \item{crown_type}{Type of the crown between paraboloid (P) and ellispoid (E). (character)}
#'        \item{h_m}{Height of the tree trunk in m. (double)}
#'        \item{hbase_m}{Crown base height of the tree (i.e. height at wich the crown start) in m. (double)}
#'        \item{hmax_m}{Height at which the crown radius is maximum. 
#'        Set it to NA in case of simple crown shapes ("P" and "E"), otherwise, your values won't be considered as hmax is automatically computed.}
#'        \item{rn_m, rs_m, re_m, rw_m}{Biggest radius of the tree crown in m (double)}
#'        \item{crown_openess}{Crown Openness of the tree (no unit) (i.e. Fraction of the energy of a light ray crossing the crown that is intercepted).
#'          Used when computing interception with a crown considered as a porous envelop. (double)}
#'        \item{crown_lad}{Leaf Area Density of the tree crown in m2/m3 (i.e. surface of leave per volume of crown, considering an homogeneous crown). 
#'          Used when computing interception with a crown considered as a turbid medium. (double)}
#'    }
#'   }
#'   \item{species_color}{named vector that contains color HEX code for each species}
#'   \item{radiation}{data.frame that contains monthly radiations.
#'    \describe{
#'      \item{month}{Unique id of the tree. (integer)}
#'      \item{Hrad}{Monthly radiation in a horizontal plane in MJ. (double)}
#'      \item{DGratio}{Ratio between monthly diffuse and global energies. (double)}
#'    }
#'   }
#'   \item{info}{A named numeric vector with site information:
#'     \describe{
#'       \item{latitude}{Latitude of the site (decimal degrees).}
#'       \item{longitude}{Longitude of the site (decimal degrees).}
#'       \item{size_x}{Horizontal extent of the site (meters).}
#'       \item{size_y}{Vertical extent of the site (meters).}
#'       \item{slope}{Mean slope of the terrain (degrees).}
#'       \item{aspect}{Aspect (degrees from north).}
#'       \item{north_to_x_cw}{Angle from north to the x-axis (clockwise, degrees).}
#'    }
#'   }
#' }
#' 
#' @source Courbaud Benoit (INRAe LESSEM Grenoble)
#' @name data_prenovel
#' @docType data
#' @usage data(data_prenovel)
#' @keywords datasets
"data_prenovel"


#' Information for IRRES1 inventory
#'
#' @description Information about IRRES1 inventory in Ardennes (Belgium) used 
#' as an example for informing the inventory zone with a core polygon.
#' 
#' @format `data_IRRES1`, a named list with 5 elements:
#'  \describe{
#'    \item{trees}{data.frame that contains for each tree, the information about 
#'      its location, species, size, size of the crown and crown information.
#'      \describe{
#'        \item{id_tree}{Unique id of the tree. (integer)}
#'        \item{species}{Species Latin name. (character)}
#'        \item{x, y}{Coordinates of the base of the tree in the forest stand in m. (double)}
#'        \item{dbh_cm}{Diameter at breast height (1.30m) of the trunk of the tree in cm. (double)}
#'        \item{crown_type}{Type of the crown between paraboloid (P) and ellispoid (E). (character)}
#'        \item{h_m}{Height of the tree trunk in m. (double)}
#'        \item{hbase_m}{Crown base height of the tree (i.e. height at wich the crown start) in m. (double)}
#'        \item{hmax_m}{Height at which the crown radius is maximum. 
#'        Set it to NA in case of simple crown shapes ("P" and "E"), otherwise, your values won't be considered as hmax is automatically computed.}
#'        \item{rn_m, rs_m, re_m, rw_m}{Biggest radius of the tree crown in m (double)}
#'        \item{crown_openess}{Crown Openness of the tree (no unit) (i.e. Fraction of the energy of a light ray crossing the crown that is intercepted).
#'          Used when computing interception with a crown considered as a porous envelop. (double)}
#'        \item{crown_lad}{Leaf Area Density of the tree crown in m2/m3 (i.e. surface of leave per volume of crown, considering an homogeneous crown). 
#'          Used when computing interception with a crown considered as a turbid medium. (double)}
#'    }
#'   }
#'   \item{species_color}{named vector that contains color HEX code for each species}
#'   \item{sensors}{data.frame that contains information about the stand sensors}
#'   \describe{
#'    \item{id_sensor}{Unique id of the sensor (integer)}
#'    \item{x, y}{Coordinates of the sensor in the stand (double)}
#'    \item{h_m}{Height of the sensor in m. (double)}
#'    \item{pacl, pacl_direct, pacl_diffuse}{Measured proportion of above canopy light on the sensor, 
#'      for respectively total, direct and diffuse energies}
#'   }
#'   \item{core_polygon}{data.frame containing vertices of the tree inventory zone, described by the coordinates (x, y) of each vertex edge.}
#'   \item{radiation}{data.frame that contains monthly radiations.
#'    \describe{
#'      \item{month}{Unique id of the tree. (integer)}
#'      \item{Hrad}{Monthly radiation in a horizontal plane in MJ. (double)}
#'      \item{DGratio}{Ratio between monthly diffuse and global energies. (double)}
#'    }
#'   }
#'   \item{info}{A named numeric vector with site information:
#'     \describe{
#'       \item{latitude}{Latitude of the site (decimal degrees).}
#'       \item{longitude}{Longitude of the site (decimal degrees).}
#'       \item{size_x}{Horizontal extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{size_y}{Vertical extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{slope}{Mean slope of the terrain (degrees).}
#'       \item{aspect}{Aspect (degrees from north).}
#'       \item{north_to_x_cw}{Angle from north to the x-axis (clockwise, degrees).}
#'    }
#'   }
#' }
#' 
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
#' @name data_IRRES1
#' @docType data
#' @usage data(data_IRRES1)
#' @keywords datasets
"data_IRRES1"


#' Information for bechefa inventory
#'
#' @description Information about bechefa marteloscope in Belgium used 
#' as an example estimating light interception using complex asymmetric crown shapes.
#' 
#' @format `data_bechefa`, a named list with 4 elements:
#'  \describe{
#'    \item{trees}{data.frame that contains for each tree, the information about 
#'      its location, species, size, size of the crown and crown information.
#'      \describe{
#'        \item{id_tree}{Unique id of the tree. (integer)}
#'        \item{species}{Species Latin name. (character)}
#'        \item{x, y}{Coordinates of the base of the tree in the forest stand in m. (double)}
#'        \item{dbh_cm}{Diameter at breast height (1.30m) of the trunk of the tree in cm. (double)}
#'        \item{crown_type}{Type of the crown between paraboloid (P) and ellispoid (E). (character)}
#'        \item{h_m}{Height of the tree trunk in m. (double)}
#'        \item{hbase_m}{Crown base height of the tree (i.e. height at wich the crown start) in m. (double)}
#'        \item{hmax_m}{Height at which the crown radius is maximum. 
#'        Set it to NA in case of simple crown shapes ("P" and "E"), otherwise, your values won't be considered as hmax is automatically computed.}
#'        \item{rn_m, rs_m, re_m, rw_m}{Biggest radius of the tree crown in m (double)}
#'        \item{crown_openess}{Crown Openness of the tree (no unit) (i.e. Fraction of the energy of a light ray crossing the crown that is intercepted).
#'          Used when computing interception with a crown considered as a porous envelop. (double)}
#'        \item{crown_lad}{Leaf Area Density of the tree crown in m2/m3 (i.e. surface of leave per volume of crown, considering an homogeneous crown). 
#'          Used when computing interception with a crown considered as a turbid medium. (double)}
#'    }
#'   }
#'   \item{species_color}{named vector that contains color HEX code for each species}
#'   \item{core_polygon}{data.frame containing vertices of the tree inventory zone, described by the coordinates (x, y) of each vertex edge.}
#'   \item{radiation}{data.frame that contains monthly radiations.
#'    \describe{
#'      \item{month}{Unique id of the tree. (integer)}
#'      \item{Hrad}{Monthly radiation in a horizontal plane in MJ. (double)}
#'      \item{DGratio}{Ratio between monthly diffuse and global energies. (double)}
#'    }
#'   }
#'   \item{info}{A named numeric vector with site information:
#'     \describe{
#'       \item{latitude}{Latitude of the site (decimal degrees).}
#'       \item{longitude}{Longitude of the site (decimal degrees).}
#'       \item{size_x}{Horizontal extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{size_y}{Vertical extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{slope}{Mean slope of the terrain (degrees).}
#'       \item{aspect}{Aspect (degrees from north).}
#'       \item{north_to_x_cw}{Angle from north to the x-axis (clockwise, degrees).}
#'    }
#'   }
#' }
#' 
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
#' @name data_bechefa
#' @docType data
#' @usage data(data_bechefa)
#' @keywords datasets
"data_bechefa"


#' Information for cloture20 inventory
#'
#' @description Information about cloture20 inventory in Belgium used 
#' as an example for a complex core polygon and demonstrating light sensors.
#' 
#' @format `data_cloture20`, a named list with 5 elements:
#'  \describe{
#'    \item{trees}{data.frame that contains for each tree, the information about 
#'      its location, species, size, size of the crown and crown information.
#'      \describe{
#'        \item{id_tree}{Unique id of the tree. (integer)}
#'        \item{species}{Species Latin name. (character)}
#'        \item{x, y}{Coordinates of the base of the tree in the forest stand in m. (double)}
#'        \item{dbh_cm}{Diameter at breast height (1.30m) of the trunk of the tree in cm. (double)}
#'        \item{crown_type}{Type of the crown between paraboloid (P) and ellispoid (E). (character)}
#'        \item{h_m}{Height of the tree trunk in m. (double)}
#'        \item{hbase_m}{Crown base height of the tree (i.e. height at wich the crown start) in m. (double)}
#'        \item{hmax_m}{Height at which the crown radius is maximum. 
#'        Set it to NA in case of simple crown shapes ("P" and "E"), otherwise, your values won't be considered as hmax is automatically computed.}
#'        \item{rn_m, rs_m, re_m, rw_m}{Biggest radius of the tree crown in m (double)}
#'        \item{crown_openess}{Crown Openness of the tree (no unit) (i.e. Fraction of the energy of a light ray crossing the crown that is intercepted).
#'          Used when computing interception with a crown considered as a porous envelop. (double)}
#'        \item{crown_lad}{Leaf Area Density of the tree crown in m2/m3 (i.e. surface of leave per volume of crown, considering an homogeneous crown). 
#'          Used when computing interception with a crown considered as a turbid medium. (double)}
#'    }
#'   }
#'   \item{species_color}{named vector that contains color HEX code for each species}
#'   \item{sensors}{data.frame that contains information about the stand sensors}
#'   \describe{
#'    \item{id_sensor}{Unique id of the sensor (integer)}
#'    \item{x, y}{Coordinates of the sensor in the stand (double)}
#'    \item{h_m}{Height of the sensor in m. (double)}
#'    \item{pacl, pacl_direct, pacl_diffuse}{Measured proportion of above canopy light on the sensor, 
#'      for respectively total, direct and diffuse energies}
#'   }
#'   \item{core_polygon}{data.frame containing vertices of the tree inventory zone, described by the coordinates (x, y) of each vertex edge.}
#'   \item{radiation}{data.frame that contains monthly radiations.
#'    \describe{
#'      \item{month}{Unique id of the tree. (integer)}
#'      \item{Hrad}{Monthly radiation in a horizontal plane in MJ. (double)}
#'      \item{DGratio}{Ratio between monthly diffuse and global energies. (double)}
#'    }
#'   }
#'   \item{info}{A named numeric vector with site information:
#'     \describe{
#'       \item{latitude}{Latitude of the site (decimal degrees).}
#'       \item{longitude}{Longitude of the site (decimal degrees).}
#'       \item{size_x}{Horizontal extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{size_y}{Vertical extent of the site (meters). Here, it is set to NA as the inventory is defined by the core polygon.}
#'       \item{slope}{Mean slope of the terrain (degrees).}
#'       \item{aspect}{Aspect (degrees from north).}
#'       \item{north_to_x_cw}{Angle from north to the x-axis (clockwise, degrees).}
#'    }
#'   }
#' }
#' 
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
#' @name data_cloture20
#' @docType data
#' @usage data(data_cloture20)
#' @keywords datasets
"data_cloture20"
