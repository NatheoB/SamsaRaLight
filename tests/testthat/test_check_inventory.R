library(testthat)

context("check_inventory")

# Load example dataset
data_prenovel <- SamsaRaLight::data_prenovel
trees <- data_prenovel$trees

test_that("valid inventory passes", {
  expect_invisible(check_inventory(trees))
})

test_that("missing required column triggers error", {
  trees_missing <- trees[, !(names(trees) %in% "x")]
  expect_error(check_inventory(trees_missing), "Missing required column")
})

test_that("non-numeric numeric columns triggers error", {
  trees_bad <- trees
  trees_bad$h_m <- as.character(trees_bad$h_m)
  expect_error(check_inventory(trees_bad), "must be numeric")
})

test_that("negative numeric values triggers error", {
  trees_bad <- trees
  trees_bad$rn_m[1] <- -1
  expect_error(check_inventory(trees_bad), "must be non-negative")
})

test_that("invalid crown_type triggers error", {
  trees_bad <- trees
  trees_bad$crown_type[1] <- "INVALID"
  expect_error(check_inventory(trees_bad), "must be one of")
})

test_that("missing hmax_m triggers error when required", {
  trees_bad <- trees
  trees_bad$crown_type[1] <- "2E"
  trees_bad$hmax_m[1] <- NA
  expect_error(check_inventory(trees_bad), "hmax_m` is mandatory for crown_type")
})

test_that("missing directional radii triggers error", {
  trees_bad <- trees
  trees_bad$crown_type[1] <- "8E"
  trees_bad$rn_m[1] <- NA
  expect_error(check_inventory(trees_bad), "All crown radii must be provided")
})

test_that("each tree must have crown_openess or crown_lad", {
  trees_bad <- trees
  trees_bad$crown_openess <- NA
  trees_bad$crown_lad <- NA
  expect_error(check_inventory(trees_bad), "must have at least one crown interception property")
})

test_that("warnings for optional columns are emitted", {
  trees_warn <- trees
  trees_warn$dbh_cm <- NA
  trees_warn$crown_openess <- NA
  trees_warn$crown_lad <- NA
  expect_warning(check_inventory(trees_warn), "Trunk interception cannot be estimated")
  expect_warning(check_inventory(trees_warn), "Porous envelope interception unavailable")
  expect_warning(check_inventory(trees_warn), "Turbid medium interception unavailable")
})

test_that("verbose = FALSE suppresses warnings/messages", {
  trees_warn <- trees
  trees_warn$dbh_cm <- NA
  trees_warn$crown_openess <- NA
  trees_warn$crown_lad <- NA
  expect_silent(check_inventory(trees_warn, verbose = FALSE))
})

# -------------------------
# Tests for symmetric crown warning
# -------------------------
test_that("symmetric crowns with differing radii trigger warning", {
  trees_sym <- trees
  # force symmetric type
  trees_sym$crown_type[1] <- "E"
  # provide differing radii
  trees_sym$rn_m[1] <- 1
  trees_sym$rs_m[1] <- 2
  trees_sym$re_m[1] <- 3
  trees_sym$rw_m[1] <- 4
  expect_warning(
    check_inventory(trees_sym),
    "Different crown radii provided for symmetric crown types"
  )
})

test_that("symmetric crowns with identical radii do not trigger warning", {
  trees_sym <- trees
  trees_sym$crown_type[1] <- "P"
  trees_sym$rn_m[1] <- 2
  trees_sym$rs_m[1] <- 2
  trees_sym$re_m[1] <- 2
  trees_sym$rw_m[1] <- 2
  expect_silent(check_inventory(trees_sym))
})

# -------------------------
# Tests for hmax_m provided but ignored
# -------------------------
test_that("hmax_m provided but ignored triggers warning for non-2E/8E crowns", {
  trees_warn <- trees
  # force non-required crown type
  trees_warn$crown_type[1] <- "E"
  trees_warn$hmax_m[1] <- 10
  expect_warning(
    check_inventory(trees_warn),
    "`hmax_m` is provided but will be ignored"
  )
})

test_that("hmax_m provided and required does not trigger ignored warning", {
  trees_ok <- trees
  # force required crown type
  trees_ok$crown_type[1] <- "2E"
  trees_ok$hmax_m[1] <- 10
  expect_silent(check_inventory(trees_ok))
})
