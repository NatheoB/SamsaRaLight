#' Example forest inventory datasets for SamsaRaLight
#'
#' @description
#' These datasets provide example forest inventories used for light interception
#' simulations with the SamsaRaLight package. Each dataset is a named list
#' with 5 elements: `trees`, `sensors`, `core_polygon`, `radiations`, and `info`.
#' 
#' @format A named list with 5 elements:
#' \describe{
#'   \item{trees}{A `data.frame` with tree-level data:
#'     \describe{
#'       \item{id_tree}{Unique id of the tree (integer).}
#'       \item{species}{Latin species name (character).}
#'       \item{x, y}{Coordinates of the base of the tree in meters (numeric).}
#'       \item{dbh_cm}{Diameter at breast height in cm (numeric).}
#'       \item{crown_type}{Crown shape type, e.g. "P", "E", or complex types like "8E" (character).}
#'       \item{h_m}{Height of the tree trunk in meters (numeric).}
#'       \item{hbase_m}{Height of crown base in meters (numeric).}
#'       \item{hmax_m}{Height at which crown radius is maximum, NA if not used (numeric).}
#'       \item{rn_m, rs_m, re_m, rw_m}{Crown radii in meters (numeric).}
#'       \item{crown_lad}{Leaf area density of the crown (m2/m3) (numeric).}
#'     }
#'   }
#'   \item{sensors}{A `data.frame` with sensor data:
#'     \describe{
#'       \item{id_sensor}{Unique sensor ID (integer).}
#'       \item{x, y}{Coordinates of the sensor in meters (numeric).}
#'       \item{h_m}{Height of the sensor in meters (numeric).}
#'       \item{pacl, pacl_direct, pacl_diffuse}{Proportion of above-canopy light measured at the sensor (numeric).}
#'     }
#'     May be `NULL` if no sensors are present.}
#'   \item{core_polygon}{A `data.frame` with vertices of the inventory polygon:
#'     \describe{
#'       \item{x, y}{Coordinates of polygon vertices in meters (numeric).}
#'     }
#'   }
#'   \item{radiations}{A `data.frame` of monthly radiation:
#'     \describe{
#'       \item{month}{Month number (integer).}
#'       \item{Hrad}{Monthly radiation in MJ (numeric).}
#'       \item{DGratio}{Diffuse-to-global radiation ratio (numeric).}
#'     }
#'   }
#'   \item{info}{A named list with site information:
#'     \describe{
#'       \item{latitude, longitude}{Coordinates of the site in decimal degrees (numeric).}
#'       \item{slope}{Mean slope in degrees (numeric).}
#'       \item{aspect}{Aspect in degrees from north (numeric).}
#'       \item{north2x}{Angle from north to x-axis clockwise in degrees (numeric).}
#'     }
#'   }
#' }
#' @usage
#' data(data_prenovel)
#' data(data_IRRES1)
#' data(data_bechefa)
#' data(data_cloture20)
#' @keywords datasets
#' @name SamsaRaLight_data
NULL

#' Prenovel inventory dataset
#'
#' @rdname SamsaRaLight_data
#' @source Courbaud Benoit (INRAe LESSEM Grenoble)
"data_prenovel"

#' IRRES1 inventory dataset
#'
#' @rdname SamsaRaLight_data
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
"data_IRRES1"

#' Bechefa inventory dataset
#'
#' @rdname SamsaRaLight_data
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
"data_bechefa"

#' Cloture20 inventory dataset
#'
#' @rdname SamsaRaLight_data
#' @source Gauthier Ligot (Gembloux Agro-Bio Tech)
"data_cloture20"
