# Install released version from CRAN
install.packages("pkgdown")

# Run once to configure your package to use pkgdown
# usethis::use_pkgdown()

# Build website
pkgdown::build_site()
