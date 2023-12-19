library(dplyr)
library(usethis)

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
    crown_type = case_match(
      species,
      c("Abies alba", "Picea abies") ~ "P",
      "Fagus sylvatica" ~ "E"
    ),
    crown_lad = 0.5,
    crown_openess = 0.2
  ) %>% 
  dplyr::select(id_tree = Id, species,
                x = X, y = Y, dbh_cm = Dbh, 
                crown_type,
                h_m = H, hbase_m = CBH, r_m = CR,
                crown_openess, crown_lad)

# Add dataset to package resources
usethis::use_data(data_trees_prenovel, overwrite = TRUE)
