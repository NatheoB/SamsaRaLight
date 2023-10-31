#' Trees dataset for Prenovel stand (Jura)
#'
#' @description Dataset containing information about Prenovel forest stand in Jura (France).
#' It contains for each trees of the plot, the information about loalisation, species, size,
#' size of the crown and crown information. All the variables needed to run SamsaRaLight.
#'
#' @format `data_trees_prenovel`
#' A data.frame with 333 rows and 12 columns:
#' \describe{
#'    \item{id_tree}{Unique id of the tree. (integer)}
#'    \item{species}{Species Latin name. (character)}
#'    \item{x, y, z}{Coordinates of the base of the tree in the forest stand in m. (double)}
#'    \item{dbh_cm}{Diameter at breast height (1.30m) of the trunk of the tree in cm. (double)}
#'    \item{height_m}{Height of the tree trunk in m. (double)}
#'    \item{cbh_m}{Crown base height of the tree (i.e. height at wich the crown start) in m. (double)}
#'    \item{cradius_m}{Biggest radius of the tree crown in m (double)}
#'    \item{crown_type}{Type of the crown between "paraboloid" (1) and "ellispoid (2). (integer)}
#'    \item{crown_lad}{Leaf Area Density of the tree crown in m2/m3 (i.e. surface of leave per volume of crown, considering an homogeneous crown). 
#'    Used when computing interception with a crown considered as a turbid medium. (double)}
#'    \item{crown_p}{Crown Openness of the tree (no unit) (i.e. Fraction of the energy of a light ray crossing the crown that is intercepted).
#'    Used when computing interception with a crown considered as a porous envelop. (double)}
#'    }
#'
#' @source Courbaud Benoit (LESSEM Grenoble)
#'
"data_trees_prenovel"