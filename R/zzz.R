# Prevent R CMD check notes for dplyr, ggplot2, sf pipelines, and internal variables
if (getRversion() >= "2.15.1") utils::globalVariables(
  c(".data",
    "x", "y", "rn_m", "rs_m", "re_m", "rw_m", "h_m", 
    "hbase_m", "hmax_m", "dbh_cm", "crown_type", "id_tree",
    "species", "crown_openness", "crown_lad",
    "X", "Y", "z", "r", "x_center", "y_center", "z_center", "id_cell",
    "month", "year", "Hrad", "DGratio",
    "added_to_fill", "view", "pos", "r_left_m", "r_right_m",
    "id_sensor", "e", "pacl", "punobs", "epot", "lci", "eunobs", "zval"
  )
)