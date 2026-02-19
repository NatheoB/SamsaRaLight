test_that("IRRES stand runs correctly", {
  
  data(data_IRRES1, package = "SamsaRaLight")
  
  inv <- data_IRRES1$trees
  sensors_xy <- create_xy_from_lonlat(data_IRRES1$sensors)$df
  info <- data_IRRES1$info
  
  expect_silent(check_coordinates(inv, verbose = FALSE))
  expect_silent(check_inventory(inv, verbose = FALSE))
  expect_silent(check_sensors(sensors_xy, verbose = FALSE))
  
  stand <- create_sl_stand(
    inv,
    cell_size = 5,
    latitude = info$latitude,
    slope = info$slope,
    aspect = info$aspect,
    north2x = info$north2x,
    sensors = sensors_xy,
    core_polygon_df = data_IRRES1$core_polygon,
    modify_polygon = "aarect",
    fill_around = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(stand, "sl_stand")
  
  radiations <- get_monthly_radiations(info$latitude, info$longitude)
  expect_silent(check_monthly_radiations(radiations, verbose = FALSE))
  
  out <- run_sl(
    stand,
    radiations,
    sensors_only = FALSE,
    detailed_output = TRUE,
    parallel_mode = FALSE, 
    verbose = FALSE
  )
  
  expect_s3_class(out, "sl_output")
  expect_true(length(out) > 0)
  
})