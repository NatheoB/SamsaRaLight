library(testthat)
library(SamsaRaLight)
library(sf)
library(dplyr)

context("create_rect_stand")

data_prenovel <- SamsaRaLight::data_prenovel
trees_inv <- data_prenovel$trees

test_that("basic valid inventory works", {
  stand <- create_rect_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    use_rect_zone = FALSE,
    fill_around = FALSE
  )
  
  expect_true(is.list(stand))
  expect_true(all(c("trees", "inv_zone_df", "inv_zone_sf", "info") %in% names(stand)))
  expect_true(is.data.frame(stand$trees))
  expect_true(all(c("x","y","dbh_cm","id_tree","added_to_fill") %in% names(stand$trees)))
  expect_true(all(stand$trees$added_to_fill == FALSE))
})

test_that("use_rect_zone = TRUE returns rectangle", {
  stand_rect <- create_rect_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    use_rect_zone = TRUE,
    fill_around = FALSE
  )
  
  expect_true(nrow(stand_rect$inv_zone_df) %in% c(4,5)) # 4 distinct corners, last might repeat
  expect_true(sf::st_geometry_type(stand_rect$inv_zone_sf)[1] == "POLYGON")
})

test_that("fill_around adds extra trees", {
  stand_filled <- create_rect_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    use_rect_zone = TRUE,
    fill_around = TRUE
  )
  
  expect_true(any(stand_filled$trees$added_to_fill))
  expect_true(nrow(stand_filled$trees) > nrow(trees_inv))
})

test_that("core polygon less than 3 vertices is ignored", {
  small_poly <- data.frame(x = c(1,2), y = c(1,2))
  
  expect_warning(
    create_rect_stand(
      trees_inv = trees_inv,
      cell_size = 5,
      core_polygon_df = small_poly,
      use_rect_zone = FALSE
    ),
    "less than 3 vertices"
  )
})

test_that("all trees remain within polygon after shifting", {
  stand <- create_rect_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    use_rect_zone = FALSE,
    fill_around = FALSE
  )
  
  coords_sf <- sf::st_as_sf(stand$trees, coords = c("x","y"))
  inside <- sf::st_intersects(coords_sf, stand$inv_zone_sf, sparse = FALSE)
  expect_true(all(inside))
})

test_that("buffering works if tree slightly outside polygon", {
  core_poly <- data.frame(
    x = range(trees_inv$x)[1:3],
    y = range(trees_inv$y)[1:3]
  )
  
  expect_warning(
    stand_buffer <- create_rect_stand(
      trees_inv = trees_inv,
      cell_size = 5,
      core_polygon_df = core_poly,
      use_rect_zone = FALSE,
      fill_around = FALSE
    ),
    "buffer"
  )
})


test_that("verbose option works", {
  expect_silent(create_rect_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    verbose = FALSE
  ))
})