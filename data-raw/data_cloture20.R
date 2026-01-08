library(dplyr)
library(usethis)

# Filepath of the SamsaraLightLoader inventory file for Cloture20 stand
inv_fp <- "data-raw/Cloture20.inv"

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
data_trees_cloture20 <- trees %>% 
  dplyr::mutate(
    species = case_match(
      SpCode,
      3 ~ "Fagus sylvatica",
      1 ~ "Quercus sp."
    ),
    crown_type = "8E") %>% 
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
data_sensors_cloture20 <- sensors %>% 
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
core_polygon_cloture20 <- core_polygon %>% 
  dplyr::select(x = X, y = Y)

  
  
  
# CREATE RADIATION DATASET ----

source("R/get_monthly_radiations.R")

# Coordinates of Cloture20 stand
longitude <- 5.20634
latitude <- 50.03617

# Fetch radiation data from PVGIS
data_rad_cloture20 <- get_monthly_radiations(latitude = latitude,
                                             longitude = longitude,
                                             start_year = 2005,
                                             end_year = 2020)


# Create colors set for species ----
species_colors_cloture20 <- c("Quercus sp." = "#CBD600",
                              "Fagus sylvatica" = "#F76045")


# FORMAT THE DATASETS AS A LIST ----
data_cloture20 <- list(
  "trees" = data_trees_cloture20,
  "species_colors" = species_colors_cloture20,
  "sensors" = data_sensors_cloture20,
  "core_polygon" = core_polygon_cloture20,
  "radiations" = data_rad_cloture20,
  "info" = list("latitude" = latitude,
                "longitude" = longitude,
                "slope" = 0,
                "aspect" = 0,
                "north2x" = 90)
)


# ADD DATASETS TO PACKAGE RESOURCES ----
usethis::use_data(data_cloture20, overwrite = TRUE)
