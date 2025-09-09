rm(list = ls())

# Libraries ----
library(dplyr)
library(microbenchmark)
library(vroom)
library(ggplot2)
library(RCapsis)
library(SamsaRaLight)

source("R/sl_run_capsis.R")

# Init for SamsaRaLight R ----
trees <- SamsaRaLight::data_trees_prenovel %>%
  dplyr::mutate(rn_m = r_m, re_m = r_m, rs_m = r_m, rw_m = r_m,
                crown_lad = case_match(species,
                                       c("Picea abies", "Abies alba") ~ 0.767,
                                       "Fagus sylvatica" ~ 0.96),
                crown_type = case_match(species,
                                        c("Picea abies", "Abies alba") ~ "P",
                                        "Fagus sylvatica" ~ "E"),
                hmax_m = NA)

monthly_rad <- data.frame(
  month = 1:12,
  Hrad = c(0,0,0,0, 533.5, 573.9, 668.9, 565.8, 418.2, 0,0,0),
  DGratio = c(0.5,0.5,0.5,0.5, 0.539146523, 0.52980456, 0.42936429, 0.435726902, 0.463048184,0.5,0.5,0.5)
)


# Compare time performance ----
microbenchmark(
  
  "r" = sl_run(
    trees, monthly_rad,
    start_day = 121, end_day = 273,
    latitude = 46, slope = 6,
    aspect = 144, north_to_x_cw = 54,
    cell_size = 10, n_cells_x = 10, n_cells_y = 10,
    turbid_medium = TRUE,
    trunk_interception = FALSE
  ),

  "samsara_rcapsis_withinit" = sl_run_samsara(capsis_folderpath = "c:/capsis",
                                       inv_fp = "C:/capsis/data/samsara2/DemoFiles/Prenovel/Prenovel01_26_Inventory_lad.txt",
                                       use_rcapsis = TRUE,
                                       java_server_started = FALSE),
  
  "lilo_rcapsis_withinit" = sl_run_lilo(capsis_folderpath = "c:/capsis4",
                                        inv_fp = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc/inv",
                                        meteo_fp = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc/meteo",
                                        export_dir = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc",
                                        run_cell = FALSE,
                                        java_server_started = FALSE),

  "samsara_cmdcapsis" = sl_run_samsara(capsis_folderpath = "c:/capsis",
                                       inv_fp = "C:/capsis/data/samsara2/DemoFiles/Prenovel/Prenovel01_26_Inventory_lad.txt",
                                       use_rcapsis = FALSE,
                                       commandfile_fp = "C:/capsis/data/samsara2/DemoFiles/Prenovel/Prenovel_CommandFile_light.txt"),
  times = 100, unit = "ms"
)

setCapsisPath("c:/capsis")
connectToCapsis()
microbenchmark(
  "samsara_rcapsis_noinit" = sl_run_samsara(capsis_folderpath = "c:/capsis",
                                            inv_fp = "C:/capsis/data/samsara2/DemoFiles/Prenovel/Prenovel01_26_Inventory_lad.txt",
                                            use_rcapsis = TRUE,
                                            java_server_started = TRUE),
  times = 100, unit = "ms"
)
shutdownClient()


setCapsisPath("c:/capsis4")
connectToCapsis()
microbenchmark(
  "lilo_rcapsis_noinit" = sl_run_lilo(capsis_folderpath = "c:/capsis4",
                                      inv_fp = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc/inv",
                                      meteo_fp = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc/meteo",
                                      export_dir = "C:/Users/nbeauchamp/Desktop/SamsaRaLight/dev/prenovel_rsc",
                                      run_cell = FALSE,
                                      java_server_started = TRUE),
  # For lilo_rcapsis and prenovel, I had set north to x angle to 90 (because disccrete values between
  # (0-90-180-270)
  times = 100, unit = "ms"
)
shutdownClient()


# Output ---- 
# (median for 100 repetitions in milliseconds)

# r                         46.4569  
# lilo_rcapsis_noinit       401.933
# samsara_rcapsis_noinit    1295.125
# samsara_cmdcapsis         2606.0432
# lilo_rcapsis_withinit     4690.8261 
# samsara_rcapsis_withinit  6753.7858 


