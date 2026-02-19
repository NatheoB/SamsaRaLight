test_that("Cloture stand runs correctly", {
  
  data(data_cloture20, package = "SamsaRaLight")
  
  inv <- data_cloture20$trees
  sensors <- data_cloture20$sensors
  info <- data_cloture20$info
  
  expect_silent(check_coordinates(inv, verbose = FALSE))
  expect_silent(check_inventory(inv, verbose = FALSE))
  expect_silent(check_sensors(sensors, verbose = FALSE))
  
  stand <- create_sl_stand(
    inv,
    cell_size = 5,
    latitude = info$latitude,
    slope = info$slope,
    aspect = info$aspect,
    north2x = info$north2x,
    sensors = sensors,
    core_polygon_df = data_cloture20$core_polygon,
    modify_polygon = "none",
    fill_around = TRUE,
    verbose = FALSE
  )
  
  expect_s3_class(stand, "sl_stand")
  
  radiations <- get_monthly_radiations(info$latitude, info$longitude)
  expect_silent(check_monthly_radiations(radiations, verbose = FALSE))
  
  out <- run_sl(
    stand,
    radiations,
    sensors_only = FALSE,
    detailed_output = FALSE,
    parallel_mode = FALSE, 
    verbose = FALSE
  )
  
  expect_s3_class(out, "sl_output")
  expect_true(length(out) > 0)
  
})