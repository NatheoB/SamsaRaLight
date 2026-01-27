library(dplyr)
library(usethis)

# CREATE TREES DATASET ----

# Filepath of Samsara inventory file for Prenovel stand
inv_fp <- "data-raw/Prenovel01_26_Inventory.txt"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)

# Get lines of tree records
tree_record_line_start <- grep("# Tree records", inv)
tree_record_line_end <- grep("# Sapling Records", inv)
trees <- inv[(tree_record_line_start+2):(tree_record_line_end-1)]

# Get and remove header
header <- trees[1]
trees <- trees[3:length(trees)]

header <- substr(header, 3, nchar(header))
header <- strsplit(header, "\t")[[1]]
header <- header[header != ""]

# Bind lines
trees <- paste(trees, collapse = "")

# Remove useless characters
trees <- gsub("\t\t",
              "", trees)

# Convert into data.frame
trees <- read.table(text=trees, col.names=header, sep="\t")

# Clean dataframe
data_trees_prenovel <- trees %>% 
  dplyr::mutate(species = case_match(
    speciesId,
    0 ~ "Picea abies",
    1 ~ "Abies alba",
    2 ~ "Fagus sylvatica"
  )) %>% 
  dplyr::mutate(
    
    # Set simple crown forms: paraboloid for conifers and ellipsoid for broadleaves
    crown_type = case_match(
      species,
      c("Abies alba", "Picea abies") ~ "P",
      "Fagus sylvatica" ~ "E"
    ),
    
    # LAD values comes from Beauchamp et al. 2025
    crown_lad = case_match(species,
                           c("Picea abies", "Abies alba") ~ 0.767,
                           "Fagus sylvatica" ~ 0.96),
    
    # Defualt crown openness value
    crown_openess = 0.2,
    
    # Here it is set to NA as we defined simple crown types ("P" and "E") for which maximum height is automatically computed
    hmax_m = NA
  ) %>% 
  dplyr::select(id_tree = Id, species,
                x = X, y = Y, dbh_cm = Dbh, 
                crown_type,
                h_m = H, hbase_m = CBH, hmax_m, 
                rn_m = CR, re_m = CR, rs_m = CR, rw_m = CR,
                crown_openness = crown_openess, crown_lad)



# CREATE RADIATION DATASET ----

source("R/get_monthly_radiations.R")

# Coordinates of Prenovel stand
longitude <- 5.82765
latitude <- 46.52666

# Fetch radiation data from PVGIS
data_rad_prenovel <- get_monthly_radiations(latitude = latitude,
                                            longitude = longitude,
                                            start_year = 2005,
                                            end_year = 2020)


# CREATE CORE POLYGON TABLE ----
core_polygon_prenovel <- data.frame(
  x = c(0, 100, 100, 0), 
  y = c(0, 0, 100, 100)
)


# FORMAT THE DATASETS AS A LIST ----
data_prenovel <- list(
  "trees" = data_trees_prenovel,
  "sensors" = NULL,
  "core_polygon" = core_polygon_prenovel,
  "radiations" = data_rad_prenovel,
  "info" = list("latitude" = latitude,
                "longitude" = longitude,
                "slope" = 6,
                "aspect" = 144,
                "north2x" = 54)
)


# ADD DATASETS TO PACKAGE RESOURCES ----
usethis::use_data(data_prenovel, overwrite = TRUE)
