library(dplyr)
library(usethis)

# Filepath of Samsara inventory file for bechefa stand
inv_fp <- "data-raw/IRRES1.inv"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)

# Get lines of tree records
tree_record_line_start <- grep("# Tree records", inv)
trees <- inv[(tree_record_line_start+1):length(inv)]

# Get and remove header
header <- trees[1]
trees <- trees[3:length(trees)]

header <- substr(header, 2, nchar(header))
header <- strsplit(header, "\t")[[1]]
header <- gsub(" ", "", header)

# Bind lines
trees <- paste(trees, collapse = "\n")

# Convert into data.frame
trees <- read.table(text=trees, col.names=header, sep="\t")

# Clean dataframe
data_trees_irres1 <- trees %>% 
  dplyr::mutate(
    species = case_match(
      SpCode,
      3 ~ "Fagus",
      41 ~ "Picea abies",
      43 ~ "Pseudotsuga menziesii",
      44 ~ "Larix",
      49 ~ "Abies alba"
    ),
    crown_type = "8E") %>% 
  dplyr::select(id_tree = Id, species,
                x = X, y = Y,
                dbh_cm = Dbh, crown_type,
                h_m = H,
                hbase_m = CBH, hmax_m = CMRH,
                rn_m = RN, rs_m = RS,
                re_m = RE, rw_m = RW,
                crown_openess = CrownOpenness, crown_lad = LAD)

# Add dataset to package resources
usethis::use_data(data_trees_irres1, overwrite = TRUE)
