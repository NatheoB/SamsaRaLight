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