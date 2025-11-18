library(dplyr)
library(usethis)

# Filepath of the SamsaraLightLoader inventory file for bechefa stand
inv_fp <- "data-raw/bechefa.inv"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)


# CREATE TREE RECORDS ----

# Get lines of tree records
tree_record_line_start <- grep("#tree", inv)
trees <- inv[(tree_record_line_start+1):length(inv)]

# Get header
header_trees <- trees[2]
header_trees <- substr(header_trees, 2, nchar(header_trees))
header_trees <- strsplit(header_trees, "\t")[[1]]
header_trees <- gsub(" ", "", header_trees)

# Remove header and bind lines from table
trees <- trees[3:length(trees)]
trees <- paste(trees, collapse = "\n")

# Convert into data.frame
trees <- read.table(text=trees, col.names=header_trees, sep="\t")

# Clean dataframe
data_trees_bechefa <- trees %>% 
  dplyr::mutate(
    species = case_match(
      species,
      2 ~ "fagus",
      3 ~ "carpinus",
      5 ~ "sorbus",
      6 ~ "picea",
      7 ~ "pseudotsuga",
      8 ~ "abies"
    ),
    crown_type = "8E") %>% 
  dplyr::select(id_tree = id, species,
                x, y,
                dbh_cm = dbh, crown_type,
                h_m = h,
                hbase_m = hbase, hmax_m = hmax,
                rn_m = rn, rs_m = rs,
                re_m = re, rw_m = rw,
                crown_openess = crownOpeness, crown_lad = lad)

  
  
# CREATE CORE POLYGON TABLE ----

# Get lines of inventory zone records
polygon_record_line_start <- grep("#inventory zone line", inv)
core_polygon <- inv[(polygon_record_line_start+1):(tree_record_line_start-2)]

# Get header
header_cp <- core_polygon[1]
header_cp <- substr(header_cp, 2, nchar(header_cp))
header_cp <- strsplit(header_cp, "\t")[[1]]
header_cp <- header_cp[1:2]

# Remove header and bind lines from table
core_polygon <- core_polygon[2:length(core_polygon)]
core_polygon <- paste(core_polygon, collapse = "\n")
core_polygon <- gsub("\t\t\t\t\t\t\t\t\t\t\t\t", "", core_polygon)

# Convert into data.frame
core_polygon <- read.table(text=core_polygon, col.names=header_cp, sep="\t")

# Clean dataframe
core_polygon_bechefa <- core_polygon

  
  
  
# CREATE RADIATION DATASET ----

source("R/get_monthly_rad.R")

# Coordinates of bechefa stand
longitude <- 5.2
latitude <- 50.04

# Fetch radiation data from PVGIS
data_rad_bechefa <- get_monthly_rad(latitude = latitude,
                                    longitude = longitude,
                                    start_year = 2005,
                                    end_year = 2020)

# Create colors set for species ----
species_colors_bechefa <- c("picea" = "#00D65C",
                            "abies" = "#45AAF7",
                            "fagus" = "#F76045",
                            "pseudotsuga" = "#BF21AB")


# FORMAT THE DATASETS AS A LIST ----
data_bechefa <- list(
  "trees" = data_trees_bechefa,
  "species_colors" = species_colors_bechefa,
  "sensors" = NULL,
  "core_polygon" = core_polygon_bechefa,
  "radiations" = data_rad_bechefa,
  "info" = c("latitude" = latitude,
             "longitude" = longitude,
             "size_x" = NA,
             "size_y" = NA,
             "slope" = 0,
             "aspect" = 0,
             "north_to_x_cw" = 90)
)


# ADD DATASETS TO PACKAGE RESOURCES ----
usethis::use_data(data_bechefa, overwrite = TRUE)
