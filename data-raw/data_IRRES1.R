library(dplyr)
library(usethis)
library(sf)

# Filepath of the SamsaraLightLoader inventory file for IRRES1 stand
inv_fp <- "data-raw/IRRES1.inv"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)


# CREATE TREE RECORDS ----

# Get lines of tree records
tree_record_line_start <- grep("# Tree records", inv)
trees <- inv[(tree_record_line_start+1):length(inv)]

# Get header
header_trees <- trees[1]
header_trees <- substr(header_trees, 2, nchar(header_trees))
header_trees <- strsplit(header_trees, "\t")[[1]]
header_trees <- gsub(" ", "", header_trees)

# Remove header and bind lines from table
trees <- trees[3:length(trees)]
trees <- paste(trees, collapse = "\n")

# Convert into data.frame
trees <- read.table(text=trees, col.names=header_trees, sep="\t")

# Clean dataframe
data_trees_IRRES1 <- trees %>% 
  dplyr::mutate(
    species = case_match(
      SpCode,
      3 ~ "Fagus sylvatica",
      41 ~ "Picea abies",
      43 ~ "Pseudtsuga menziesii",
      44 ~ "Larix sp.",
      49 ~ "Abies alba"
    ),
    crown_type = case_when(
      species == "Fagus sylvatica" ~ "E",
      TRUE ~ "P"
    ),
    CMRH = NA) %>% 
  dplyr::select(id_tree = Id, species,
                x = X, y = Y,
                dbh_cm = Dbh, crown_type,
                h_m = H,
                hbase_m = CBH, hmax_m = CMRH,
                rn_m = RN, rs_m = RS,
                re_m = RE, rw_m = RW,
                crown_openness = CrownOpenness, crown_lad = LAD)


# GET SENSORS RECORD TABLE ----

# Get lines of sensors records
sensors_record_line_start <- grep("# Sensor records", inv)
sensors <- inv[(sensors_record_line_start+1):(tree_record_line_start-3)]

# Get header
header_sensors <- sensors[1]
header_sensors <- substr(header_sensors, 2, nchar(header_sensors))
header_sensors <- strsplit(header_sensors, "\t")[[1]]
header_sensors <- gsub(" ", "", header_sensors)

# Remove header and bind lines from table
sensors <- sensors[2:length(sensors)]
sensors <- paste(sensors, collapse = "\n")

# Convert into data.frame
sensors <- read.table(text=sensors, col.names=header_sensors, sep="\t")

# Clean dataframe
data_sensors_IRRES1 <- sensors %>% 
  dplyr::select(
    id_sensor = ID,
    x, y, h_m = z,
    pacl = PACLtotal,
    pacl_direct = PACLdirect,
    pacl_diffuse = PACLdiffuse
  )


  
  
# CREATE CORE POLYGON TABLE ----

# Get lines of inventory zone records
polygon_record_line_start <- grep("# Inventory zone records", inv)
core_polygon <- inv[(polygon_record_line_start+1):(sensors_record_line_start-3)]

# Get header
header_cp <- core_polygon[1]
header_cp <- substr(header_cp, 3, nchar(header_cp))
header_cp <- strsplit(header_cp, "\t")[[1]]
header_cp <- gsub(" ", "", header_cp)

# Remove header and bind lines from table
core_polygon <- core_polygon[2:length(core_polygon)]
core_polygon <- paste(core_polygon, collapse = "\n")

# Convert into data.frame
core_polygon <- read.table(text=core_polygon, col.names=header_cp, sep="\t")

# Clean dataframe
core_polygon_IRRES1 <- core_polygon %>% 
  dplyr::select(x = X, y = Y)


# Convert all coordinates into lat/long ----
# FOR THE SAKE OF THE TUTORIAL 3
# Initially, coordinates are in Belgian Lambert 72 (EPSG:31370)

## Tree inventory ----
trees_sf <- st_as_sf(data_trees_IRRES1,
                     coords = c("x", "y"),
                     crs = 31370)   # Belgian Lambert 72

trees_utm <- st_transform(trees_sf, 4326) # WGS84

trees_lonlat <- st_coordinates(trees_utm)
data_trees_IRRES1$lon <- trees_lonlat[,1]
data_trees_IRRES1$lat <- trees_lonlat[,2]

# Here, when I convert it back to Lambert 72, we retrieve the same initial x and y coordinates 
# trees_sf <- st_as_sf(data_trees_IRRES1,
#                      coords = c("lon", "lat"),
#                      crs = 4326)   # WGS84
# 
# trees_utm <- st_transform(trees_sf, 31370)
# 
# xy <- st_coordinates(trees_utm)
# data_trees_IRRES1$xnew <- xy[,1]
# data_trees_IRRES1$ynew <- xy[,2]

data_trees_IRRES1 <- subset(data_trees_IRRES1, select = -c(x, y))



## Sensors ----
sensors_sf <- st_as_sf(data_sensors_IRRES1,
                       coords = c("x", "y"),
                       crs = 31370)   # Belgian Lambert 72

sensors_utm <- st_transform(sensors_sf, 4326) # WGS84

sensors_lonlat <- st_coordinates(sensors_utm)
data_sensors_IRRES1$lon <- sensors_lonlat[,1]
data_sensors_IRRES1$lat <- sensors_lonlat[,2]

data_sensors_IRRES1 <- subset(data_sensors_IRRES1, select = -c(x, y))



## Core polygon ----
polygon_sf <- st_as_sf(core_polygon_IRRES1,
                       coords = c("x", "y"),
                       crs = 31370)   # Belgian Lambert 72

polygon_utm <- st_transform(polygon_sf, 4326) # WGS84

polygon_lonlat <- st_coordinates(polygon_utm)
core_polygon_IRRES1$lon <- polygon_lonlat[,1]
core_polygon_IRRES1$lat <- polygon_lonlat[,2]

core_polygon_IRRES1 <- subset(core_polygon_IRRES1, select = -c(x, y))

  
  
# CREATE RADIATION DATASET ----

source("R/get_monthly_radiations.R")

# Coordinates of IRRES1 stand
longitude <- 5.957553
latitude <- 50.32875

# Fetch radiation data from PVGIS
data_rad_IRRES1 <- get_monthly_radiations(latitude = latitude,
                                          longitude = longitude,
                                          start_year = 2005,
                                          end_year = 2020)



# FORMAT THE DATASETS AS A LIST ----
# For the sake of the tutorial 3, set the core polygon to NULL
data_IRRES1 <- list(
  "trees" = data_trees_IRRES1,
  "sensors" = data_sensors_IRRES1,
  "core_polygon" = NULL,
  "radiations" = data_rad_IRRES1,
  "info" = list("latitude" = latitude,
                "longitude" = longitude,
                "slope" = 6.489,
                "aspect" = 43.536,
                "north2x" = 90)
)


# ADD DATASETS TO PACKAGE RESOURCES ----
usethis::use_data(data_IRRES1, overwrite = TRUE)
