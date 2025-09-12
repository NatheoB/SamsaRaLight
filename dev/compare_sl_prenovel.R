#### CAPSIS ENERGY: SAMRARA ####

# install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type="source")
# install.packages("https://sourceforge.net/projects/rcapsis.capsisbridge.p/files/latest", repos = NULL,  type="source")

library(RCapsis)
library(dplyr)
library(ggplot2)
library(data.table)
library(SamsaRaLight)

source("R/sl_run_capsis.R")

out_capsis <- sl_run_samsara("c:/capsis4",
                             "C:/capsis4/data/samsara2/DemoFiles/Prenovel/Prenovel01_26_Inventory.txt",
                             use_rcapsis = TRUE)


#### R ENERGY ####
trees <- SamsaRaLight::data_trees_prenovel %>%
                         dplyr::mutate(rn_m = r_m,
                                       re_m = r_m,
                                       rs_m = r_m,
                                       rw_m = r_m,
                                       crown_lad = case_match(species,
                                                               c("Picea abies", "Abies alba") ~ 0.767,
                                                               "Fagus sylvatica" ~ 0.96),
                                       crown_type = case_match(species,
                                                               c("Picea abies", "Abies alba") ~ "P",
                                                               "Fagus sylvatica" ~ "E"),
                                       hmax_m = case_match(crown_type,
                                                           "E" ~ hbase_m + 1/2*(h_m - hbase_m),
                                                           "P" ~ hbase_m)) %>%
                         dplyr::select(-r_m)

monthly_rad <- data.frame(
  month = 1:12,
  Hrad = c(0,0,0,0, 533.5, 573.9, 668.9, 565.8, 418.2, 0,0,0),
  DGratio = c(0.5,0.5,0.5,0.5, 0.539146523, 0.52980456, 0.42936429, 0.435726902, 0.463048184,0.5,0.5,0.5)
)

out_r <- sl_run(
  trees, monthly_rad,
  start_day = 121, end_day = 273,
  latitude = 46, slope = 6, 
  aspect = 144, north_to_x_cw = 54,
  cell_size = 10, n_cells_x = 10, n_cells_y = 10,
  turbid_medium = TRUE,
  trunk_interception = FALSE
)

microbenchmark::microbenchmark(
  "sl" = sl_run(
    trees, monthly_rad,
    start_day = 121, end_day = 273,
    latitude = 46, slope = 6, 
    aspect = 144, north_to_x_cw = 54,
    cell_size = 10, n_cells_x = 10, n_cells_y = 10,
    turbid_medium = TRUE,
    trunk_interception = FALSE
  ))


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

compare_light <- dplyr::left_join(capsis_light, r_light, by = "id_tree", suffix = c(".capsis", ".r")) %>%
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

