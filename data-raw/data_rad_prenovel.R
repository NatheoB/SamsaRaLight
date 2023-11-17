library(dplyr)
library(usethis)
source("R/get_monthly_rad.R")

# Coordinates of Prenovel stand
longitude <- 5.82765
latitude <- 46.52666

# Fetch radiation data from PVGIS
data_rad_prenovel <- get_monthly_rad(latitude = latitude,
                                     longitude = longitude,
                                     start_year = 2005,
                                     end_year = 2020)

# Add dataset to package resources
usethis::use_data(data_rad_prenovel, overwrite = TRUE)