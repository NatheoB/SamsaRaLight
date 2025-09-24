library(dplyr)
library(usethis)

# Filepath of Samsara inventory file for bechefa stand
inv_fp <- "data-raw/IRRES1.inv"

# Load inventory file
con <- file(inv_fp, "r")
inv <- readLines(con)
close(con)

# Get lines from sensor records but before tree records
sensor_record_line_start <- grep("# Sensor records", inv)
tree_record_line_start <- grep("# Tree records", inv)
sensors <- inv[(sensor_record_line_start+1):tree_record_line_start]

# Get and remove header
header <- sensors[1]
sensors <- sensors[2:length(sensors)]

header <- substr(header, 2, nchar(header))
header <- strsplit(header, "\t")[[1]]
header <- gsub(" ", "", header)

# Bind lines
sensors <- paste(sensors, collapse = "\n")

# Convert into data.frame
sensors <- read.table(text=sensors, col.names=header, sep="\t")

# Clean dataframe
data_sensors_irres1 <- sensors %>% 
  dplyr::select(id_sensor = ID,
                x, y,
                h_m = z,
                pacl = PACLtotal,
                pacl_direct = PACLdirect,
                pacl_diffuse = PACLtotal)

# Add dataset to package resources
usethis::use_data(data_sensors_irres1, overwrite = TRUE)
