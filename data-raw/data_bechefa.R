library(dplyr)
library(vroom)
library(usethis)

# CREATE TREE RECORDS ----
# Here, we rotated the bechefa stand to store it rotated for the sake of the third tutorial
data_trees_bechefa <- vroom::vroom("data-raw/bechefa_trees.csv") %>% 
  dplyr::select(-added_to_fill) %>% 
  dplyr::mutate(crown_type = "8E")

  
# CREATE CORE POLYGON TABLE ----
core_polygon_bechefa <- data.frame(
  x = c(0, 185, 185, 0),
  y = c(80, 80, 0, 0)
)

  
  
  
# CREATE RADIATION DATASET ----

source("R/get_monthly_radiations.R")

# Coordinates of bechefa stand
longitude <- 5.2
latitude <- 50.04

# Fetch radiation data from PVGIS
data_rad_bechefa <- get_monthly_radiations(latitude = latitude,
                                           longitude = longitude,
                                           start_year = 2005,
                                           end_year = 2020)


# FORMAT THE DATASETS AS A LIST ----
data_bechefa <- list(
  "trees" = data_trees_bechefa,
  "sensors" = NULL,
  "core_polygon" = core_polygon_bechefa,
  "radiations" = data_rad_bechefa,
  "info" = list("latitude" = latitude,
                "longitude" = longitude,
                "slope" = 0,
                "aspect" = 0,
                "north2x" = 68)
)


# ADD DATASETS TO PACKAGE RESOURCES ----
usethis::use_data(data_bechefa, overwrite = TRUE)
