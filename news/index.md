# Changelog

## SamsaRaLight 1.1.0

- Crown interceptions aligned with the Capsis SamsaraLight library:
  ray/crown part interceptions with an odd number of interception points
  are ignored for all crown shapes (previously undefined behaviour for
  paraboloids), with the Capsis rounding tolerance for crown fractions
  (`8E`, `4P`) and the Capsis particular cases.
- Unexpected interceptions are no longer printed in the console from the
  C++ code (not thread-safe); they are counted and reported in a single
  R warning.
- Trunks: correct position (previously uninitialised) and trunk cylinder
  from the ground to the crown base, as in Capsis (previously up to the
  tree top). Results change when `trunk_interception = TRUE` (default in
  [`run_sl()`](https://natheob.github.io/SamsaRaLight/reference/run_sl.md)).
- OpenMP parallelisation is now enabled on Linux and macOS (new
  `src/Makevars`).
- Trees or sensors exactly on the upper limits of the stand are assigned
  to the right grid cell.
- Symmetric crowns (`P`, `E`, `2E`): the crown radius is now the mean of
  the provided radii, as announced by
  [`check_inventory()`](https://natheob.github.io/SamsaRaLight/reference/check_inventory.md)
  (previously only `rw_m` was used). Results change for these crowns
  when their radii differ.
- Removed the link to the Shiny application from the website.
- Documentation fixes.

## SamsaRaLight 1.0.2

- New tree output variable `rci` (ray competition index): the proportion
  of the rays intercepted by a tree crown that were already intercepted
  by another crown (`rci_direct` and `rci_diffuse` in the detailed
  output).
- [`plot.sl_output()`](https://natheob.github.io/SamsaRaLight/reference/plot.sl_output.md):
  new option `what_trees = "compet_rays"` to map the RCI.
- Output object reorganised: simulation parameters are now in a
  top-level `params` element (including `detailed_output`), and `input`
  only contains the `sl_stand` and `monthly_radiations` objects.
- [`validate_sl_output()`](https://natheob.github.io/SamsaRaLight/reference/validate_sl_output.md)
  now checks all main output columns (`epot`, `lci`, `eunobs`, `rci`,
  `pacl`, `punobs`).
- Installation instructions now use
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html).
- Updated vignettes and documentation.

## SamsaRaLight 1.0.1

- Forced Rcpp version \>= 1.1.0 for the package (needed for compiling on
  r-devel windows and mac)
