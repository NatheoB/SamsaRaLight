#### CAPSIS ENERGY ####

# install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type="source")
# install.packages("https://sourceforge.net/projects/rcapsis.capsisbridge.p/files/latest", repos = NULL,  type="source")

library(RCapsis)
library(dplyr)
library(ggplot2)
library(data.table)
library(SamsaRaLight)


# LILO ----

source("R/sl_run_capsis.R")

out_capsis <- sl_run_lilo(capsis_folderpath = "c:/capsis4",
                          inv_fp = file.path(getwd(), "../dev/bechefa_rsc/inv"), 
                          meteo_fp = file.path(getwd(), "../dev/bechefa_rsc/meteo"), 
                          export_dir = file.path(getwd(), "../dev/bechefa_rsc"))


# R ----
rad <- data.frame(
  month = 1:12,
  Hrad = c(71.2776, 127.8588, 312.6255, 415.0542, 533.9946, 421.6902,
           425.898, 345.9888, 302.904, 181.0038, 103.2462, 64.9182),
  DGratio = 0.5)

out_r <- sl_run(data_trees_bechefa, rad,
                latitude = 50.04,
                start_day = 91, end_day = 304,
                cell_size = 5, 
                n_cells_x = 40, n_cells_y = 40,
                trunk_interception = FALSE,
                soc = FALSE,
                height_anglemin = 10,
                direct_startoffset = 7.5,
                direct_anglestep = 15,
                diffuse_anglestep = 15)

microbenchmark::microbenchmark(
  "sl" = sl_run(data_trees_bechefa, rad,
                latitude = 50.04,
                start_day = 91, end_day = 304,
                cell_size = 5, 
                n_cells_x = 40, n_cells_y = 40,
                trunk_interception = FALSE,
                soc = FALSE,
                height_anglemin = 10,
                direct_startoffset = 7.5,
                direct_anglestep = 15,
                diffuse_anglestep = 15))


#### COMPARE CAPSIS AND R ####
capsis_light <- out_capsis$trees %>%
  dplyr::arrange(desc(e)) %>%
  dplyr::mutate(rank_e = row_number()) %>%
  dplyr::arrange(desc(epot)) %>%
  dplyr::mutate(rank_epot = row_number())

r_light <- out_r$trees %>%
  as_tibble() %>%
  dplyr::select(id_tree, epot, e) %>%
  dplyr::arrange(desc(e)) %>%
  dplyr::mutate(rank_e = row_number()) %>%
  dplyr::arrange(desc(epot)) %>%
  dplyr::mutate(rank_epot = row_number())

compare_light <- dplyr::full_join(capsis_light, r_light, 
                                  by = "id_tree", 
                                  suffix = c(".capsis", ".r")) %>%
  dplyr::mutate(diff_epot = epot.capsis - epot.r)


# COMPARE POTENTIAL LIGHT
ggplot(compare_light, aes(y = rank_epot.r, x = rank_epot.capsis)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")

ggplot(compare_light, aes(y = epot.r, x = epot.capsis)) +
  geom_point() +
  geom_smooth(color = 'red', method = "lm")


# COMPARE LIGHT
ggplot(compare_light, aes(y = rank_e.r, x = rank_e.capsis)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")

ggplot(compare_light, aes(y = e.r, x = e.capsis)) +
  geom_point() +
  geom_smooth(color = 'red', method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "blue")

# all equals ?
all.equal(r_light, capsis_light)


# COMPARE CELLS
compare_cells <- dplyr::left_join(out_capsis$cells, out_r$cells,
                                  by = c("x_center", "y_center"),
                                  suffix = c(".capsis", ".r"))

ggplot(compare_cells, aes(y = e.r, x = e.capsis)) +
  geom_point() +
  geom_smooth(color = 'red', method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "blue")

ggplot(compare_cells, aes(y = erel.r, x = erel.capsis)) +
  geom_point() +
  geom_smooth(color = 'red', method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "blue")




