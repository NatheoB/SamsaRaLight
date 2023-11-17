library(devtools)
library(usethis)

# Run all tests
remove.packages("SamsaRaLight")
devtools::clean_dll()
.rs.restartR()


rm(list = ls())

devtools::load_all('.')
devtools::document('.')
devtools::test()
devtools::check()
