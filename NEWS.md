# SamsaRaLight (development version)

* New functions `precompute_sl_geometry()` and `estimate_pacl_from_geometry()` to estimate the light on sensors and cells many times for the same stand with different crown properties (e.g. calibration of crown LAD). The part of the ray tracing that does not depend on crown LAD or openness (crowns crossed by each ray, path lengths, trunks) is computed once; the light on sensors and cells is then estimated for any crown LAD or openness, extinction coefficient, clumping factor and trunk interception, with the same results as `run_sl_advanced()` up to rounding errors and 100 to 500 times faster on the package datasets. Light intercepted by trees cannot be estimated from the geometry.
* Internal: argument checks, OpenMP set-up and the warning on unexpected interceptions are shared by `run_sl_advanced()` and `precompute_sl_geometry()`; the C++ function `sl_run_rcpp()` is renamed `run_sl_rcpp()`.
* Documentation: `run_sl()` inherits its parameters from `run_sl_advanced()`, whose parameter descriptions were harmonised.
* Energy of sensors on sloping stands aligned with Capsis (`SLLighterLegacy`): the energy of a sensor or cell is now the sum of the energies reaching it along the traced rays. It was previously the energy above the canopy minus the intercepted energy, which counted as received the energy of the sky directions that are not traced because they are hidden by the slope (rays are traced only if they reach the slope with an angle above `height_anglemin`), and missed the low rays traced below this angle on the downhill side. Results change only for sensors (horizontal plane) on sloping stands, by a constant per stand (e.g. PACL −0.0003 for `data_IRRES1`, up to −0.12 on steep stands in our tests); cells, trees and flat stands are unchanged. As in Capsis, relative energies are capped at 1 (`pacl`, `pacl_direct` and `pacl_diffuse` for sensors; `pacl_direct` and `pacl_diffuse` for cells). Same change in `estimate_pacl_from_geometry()`.
* The geometry computed by `precompute_sl_geometry()` is stored as sparse matrices (new dependency on the package `Matrix`): path lengths of the rays across the crowns, crowns crossed, and energies of the rays reaching each target, with the rays intercepted by a trunk kept apart. `estimate_pacl_from_geometry()` computes the light on sensors and cells from these matrices in R (sparse matrix products), which is fast enough to be called many times, e.g. on many stands in a calibration. Geometries are not compatible between development versions: compute them again with `precompute_sl_geometry()`.

# SamsaRaLight 1.1.0

Corrections following a review of the package code by Claude Code, in comparison with the Capsis SamsaraLight Java library.

* Crown interceptions aligned with the Capsis SamsaraLight library: ray/crown part interceptions with an odd number of interception points are ignored for all crown shapes (previously undefined behaviour for paraboloids), with the Capsis rounding tolerance for crown fractions (`8E`, `4P`) and the Capsis particular cases.
* Unexpected interceptions are no longer printed in the console from the C++ code (not thread-safe); they are counted and reported in a single R warning.
* Trunks: correct position (previously uninitialised) and trunk cylinder from the ground to the crown base, as in Capsis (previously up to the tree top). Results change when `trunk_interception = TRUE` (default in `run_sl()`).
* Trunk interception in the energy balance aligned with Capsis (`SLLighterLegacy`): a trunk intercepts all the remaining energy of the ray, but the crown parts behind it are still counted in the potential energy of their tree (they were previously skipped), and the trunk no longer reduces the potential energy of its own crown. With porous envelopes, a trunk behind a part of its own crown now also stops the ray (previously it let the ray through). Interceptions at equal distance keep their order (stable sort, as in Capsis). This changes `epot`, `lci` and `rci` when `trunk_interception = TRUE`.
* OpenMP parallelisation is now enabled on Linux and macOS (new `src/Makevars`).
* Crown fractions (`8E`, `4P`): paths are now computed exactly when the target lies on a cut plane of the fraction (e.g. trees located exactly on a cell-centre coordinate). A ray running along a cut plane is counted once instead of twice, a target within the rounding tolerance of a cut plane is placed on it, and a cut-plane point coinciding with the target is no longer counted twice (which dropped the interception). Departure from Capsis, which has the same rounding tolerance but not these corrections.
* Trees or sensors exactly on the upper limits of the stand are assigned to the right grid cell.
* Symmetric crowns (`P`, `E`, `2E`): the crown radius is now the mean of the provided radii, as announced by `check_inventory()` (previously only `rw_m` was used). Results change for these crowns when their radii differ.
* `crown_lad` is no longer required in the tree inventory when crowns are porous envelopes (`turbid_medium = FALSE`).
* Removed a deprecated OpenMP construct (`master`) reported as a warning by recent compilers.
* Removed the link to the Shiny application from the website.
* Documentation fixes.

# SamsaRaLight 1.0.2

* New tree output variable `rci` (ray competition index): the proportion of the rays intercepted by a tree crown that were already intercepted by another crown (`rci_direct` and `rci_diffuse` in the detailed output).
* `plot.sl_output()`: new option `what_trees = "compet_rays"` to map the RCI.
* Output object reorganised: simulation parameters are now in a top-level `params` element (including `detailed_output`), and `input` only contains the `sl_stand` and `monthly_radiations` objects.
* `validate_sl_output()` now checks all main output columns (`epot`, `lci`, `eunobs`, `rci`, `pacl`, `punobs`).
* Installation instructions now use `pak::pkg_install()`.
* Updated vignettes and documentation.

# SamsaRaLight 1.0.1

* Forced Rcpp version >= 1.1.0 for the package (needed for compiling on r-devel windows and mac)