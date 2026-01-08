library(testthat)
library(SamsaRaLight)
library(sf)

context("prepare_input_stand")

data_prenovel <- SamsaRaLight::data_prenovel
trees_inv <- data_prenovel$trees

test_that("prepare_input_stand works with verbose = FALSE", {
  stand <- prepare_input_stand(
    trees_inv = trees_inv,
    cell_size = 5,
    verbose = FALSE
  )
  
  expect_true(is.list(stand))
  expect_true(all(c("trees","inv_zone_df","inv_zone_sf","info") %in% names(stand)))
})

test_that("prepare_input_stand shows messages with verbose = TRUE", {
  expect_message(
    prepare_input_stand(
      trees_inv = trees_inv,
      cell_size = 5,
      verbose = TRUE
    )
  )
})
