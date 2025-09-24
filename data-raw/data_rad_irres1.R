library(dplyr)
library(usethis)
source("R/get_monthly_rad.R")

# Coordinates of IRRES1 stand (same as bechefa)
longitude <- 5.2
latitude <- 50.04

# Fetch radiation data from PVGIS
data_rad_irres1 <- get_monthly_rad(latitude = latitude,
                                    longitude = longitude,
                                    start_year = 2005,
                                    end_year = 2020)

# Add dataset to package resources
usethis::use_data(data_rad_irres1, overwrite = TRUE)
