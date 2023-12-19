library(dplyr)
library(usethis)

# Filepath of Samsara inventory file for bechefa stand
inv_fp <- "data-raw/bechefa.inv"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)

# Get lines of tree records
tree_record_line_start <- grep("#tree", inv)
trees <- inv[(tree_record_line_start+2):length(inv)]

# Get and remove header
header <- trees[1]
trees <- trees[2:length(trees)]

header <- substr(header, 2, nchar(header))
header <- strsplit(header, "\t")[[1]]
header <- header[header != ""]

# Bind lines
trees <- paste(trees, collapse = "\n")

# Convert into data.frame
trees <- read.table(text=trees, col.names=header, sep="\t")

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

# Add dataset to package resources
usethis::use_data(data_trees_bechefa, overwrite = TRUE)
