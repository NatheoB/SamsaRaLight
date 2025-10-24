library(devtools)
library(usethis)
library(pkgdown)

# Run all tests
remove.packages("SamsaRaLight")
devtools::clean_dll()
.rs.restartR()

rm(list = ls())

devtools::load_all('.')
devtools::document('.')
devtools::test()
devtools::check()

pkgdown::build_site()
usethis::use_pkgdown_github_pages()
