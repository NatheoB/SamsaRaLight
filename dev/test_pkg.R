library(devtools)
library(usethis)

# Run all tests
rm(list = ls())
remove.packages("SamsaRaLight")

devtools::load_all('.')
devtools::document('.')
devtools::test()
devtools::check()
