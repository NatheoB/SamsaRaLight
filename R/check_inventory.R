#' Check the format and validity of a tree inventory data.frame
#'
#' This function checks whether a tree inventory data.frame is correctly
#' formatted to be used as input for the ray-tracing model. It verifies
#' the presence, type, and validity of mandatory and optional variables
#' describing the geometry and attributes of trees within a forest stand.
#'
#' @param trees_inv A data.frame with one row per tree and the following columns:
#' \describe{
#'   \item{id_tree}{Unique identifier of the tree (numeric or character, no duplicates).}
#'   \item{x}{X position of the tree within the stand (numeric, meters, planar coordinates).
#'     Optional if \code{lon} and \code{lat} are provided.}
#'   \item{y}{Y position of the tree within the stand (numeric, meters, planar coordinates).
#'     Optional if \code{lon} and \code{lat} are provided.}
#'   \item{lon}{Longitude of the tree location (numeric, decimal degrees).
#'     Optional; if provided with \code{lat}, coordinates are converted internally
#'     to a UTM planar coordinate system in \code{create_sl_stand()}.}
#'   \item{lat}{Latitude of the tree location (numeric, decimal degrees).
#'     Optional; see \code{lon}.}
#'   \item{species}{Species name (character).}
#'   \item{dbh_cm}{Diameter at breast height (1.3 m, in cm).}
#'   \item{crown_type}{Type of crown geometry. One of \code{"E"}, \code{"P"},
#'     \code{"2E"}, \code{"8E"}, or \code{"4P"}.}
#'   \item{h_m}{Total tree height (numeric, meters).}
#'   \item{hbase_m}{Height of the crown base (numeric, meters).}
#'   \item{hmax_m}{Height of the maximum crown radius (numeric, meters).
#'     Required only if at least one tree has a crown type \code{"2E"} or \code{"8E"}.
#'     For other crown types, the column is optional and the value is internally computed.}
#'   \item{rn_m}{Crown radius toward North (numeric, meters).}
#'   \item{rs_m}{Crown radius toward South (numeric, meters).}
#'   \item{re_m}{Crown radius toward East (numeric, meters).}
#'   \item{rw_m}{Crown radius toward West (numeric, meters).}
#'   \item{crown_openness}{Crown openness (unitless), optional if turbid medium interception.}
#'   \item{crown_lad}{Leaf Area Density (m² m⁻³), optional if porous envelope interception.}
#' }
#'
#' @details
#' The function performs the following checks and validations:
#' \describe{
#'   \item{1}{Ensures \code{trees_inv} is a non-empty data.frame with all required columns.}
#'   \item{2}{Checks that \code{id_tree} values are unique.}
#'   \item{3}{Checks that either planar coordinates (\code{x}, \code{y}) or geographic coordinates
#'            (\code{lon}, \code{lat}) are provided. If only \code{lon} and \code{lat} are supplied, they are converted to planar UTM coordinates in \code{create_sl_stand()}.
#'            The \code{plot_inventory()} function requires planar coordinates and cannot be used before this conversion.}
#'   \item{4}{Validates numeric columns (\code{dbh_cm}, \code{h_m}, \code{hbase_m}, \code{rn_m}, \code{rs_m}, \code{re_m}, \code{rw_m}) are numeric and non-negative.}
#'   \item{5}{Verifies that \code{crown_type} values are one of \code{"E"}, \code{"P"}, \code{"2E"}, \code{"8E"}, or \code{"4P"}.}
#'   \item{6}{Ensures crown radii are present according to crown type.}
#'   \item{7}{Checks that \code{hmax_m} is provided when required and lies between \code{hbase_m} and \code{h_m}.}
#'   \item{8}{Ensures \code{hbase_m < h_m}.}
#'   \item{9}{Verifies that each tree has at least one crown interception property defined.}
#'   \item{10}{Provides informative error messages and warnings for all invalid conditions.}
#' }
#'
#' @param verbose Logical; if \code{TRUE}, informative messages and warnings are printed.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#'
#' @examples
#' \dontrun{
#' # Using the example dataset from the package
#' data_prenovel <- SamsaRaLight::data_prenovel
#' trees <- data_prenovel$trees
#'
#' # Check the inventory
#' check_inventory(trees)
#'
#' # Quiet mode
#' check_inventory(trees, verbose = FALSE)
#' }
#'
#' @export
check_inventory <- function(trees_inv, verbose = TRUE) {
  
  ## ---- basic structure ------------------------------------------------------
  if (!is.data.frame(trees_inv)) stop("`trees_inv` must be a data.frame.", call. = FALSE)
  if (nrow(trees_inv) == 0) stop("`trees_inv` must contain at least one tree.", call. = FALSE)
  
  ## ---- check coordinates ----------------------------------------------------
  needs_conversion <- check_coordinates(trees_inv, verbose = F)
  
  if (needs_conversion) {
    message(
      "Geographic coordinates (`lon`, `lat`) detected. They will be converted ",
      "to planar UTM coordinates in `create_sl_stand()`. ",
      "If you want to use `plot_inventory()` before that, ",
      "convert the coordinates yourself using `convert_xy_from_lonlat()`."
    )
  }
  
  ## ---- required columns -----------------------------------------------------
  required_cols <- c("id_tree",
                     "species", "crown_type",
                     "dbh_cm", "h_m", "hbase_m", 
                     "rn_m", "rs_m", "re_m", "rw_m")
  missing_cols <- setdiff(required_cols, names(trees_inv))
  if (length(missing_cols) > 0) stop(
    "Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE
  )
  
  ## ---- id uniqueness --------------------------------------------------------
  if (anyDuplicated(trees_inv$id_tree)) stop(
    "`id_tree` must contain unique values (no duplicates).", call. = FALSE
  )
  
  ## ---- numeric checks -------------------------------------------------------
  numeric_cols <- c("dbh_cm", "h_m", "hbase_m", "rn_m", "rs_m", "re_m", "rw_m")
  non_numeric <- numeric_cols[!vapply(trees_inv[numeric_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) stop(
    "The following columns must be numeric: ", paste(non_numeric, collapse = ", "), call. = FALSE
  )
  
  numeric_pos_cols <- c("dbh_cm", "h_m", "hbase_m", "rn_m", "rs_m", "re_m", "rw_m")
  if (any(trees_inv[numeric_pos_cols] < 0, na.rm = TRUE)) stop(
    "The following columns must be non-negative: ", paste(non_numeric, collapse = ", "), call. = FALSE
  )
  
  ## ---- crown type -----------------------------------------------------------
  allowed_crown_types <- c("E", "P", "2E", "8E", "4P")
  if (!all(trees_inv$crown_type %in% allowed_crown_types)) stop(
    "`crown_type` must be one of: ", paste(allowed_crown_types, collapse = ", "), call. = FALSE
  )
  
  ## ---- conditional crown radii ---------------------------------------------
  symmetric_types <- c("E", "P", "2E")
  directional_types <- c("8E", "4P")
  radii_cols <- c("rn_m", "rs_m", "re_m", "rw_m")
  need_radii <- any(trees_inv$crown_type %in% c(symmetric_types, directional_types))
  if (need_radii) {
    missing_radii_cols <- setdiff(radii_cols, names(trees_inv))
    if (length(missing_radii_cols) > 0) stop(
      "The following crown radius column(s) are required: ", paste(missing_radii_cols, collapse = ", "), call. = FALSE
    )
  }
  
  if (all(radii_cols %in% names(trees_inv))) {
    if (!all(vapply(trees_inv[radii_cols], is.numeric, logical(1)))) stop(
      "Crown radii (rn_m, rs_m, re_m, rw_m) must be numeric.", call. = FALSE
    )
    if (any(trees_inv[radii_cols] < 0, na.rm = TRUE)) stop(
      "Crown radii must be non-negative.", call. = FALSE
    )
  }
  
  # symmetric crowns
  symmetric_rows <- trees_inv$crown_type %in% symmetric_types
  if (any(symmetric_rows)) {
    radii_mat <- as.matrix(trees_inv[symmetric_rows, radii_cols, drop = FALSE])
    different_radii <- apply(radii_mat, 1, function(x) length(unique(x[!is.na(x)])) > 1)
    if (any(different_radii) && verbose) warning(
      "Different crown radii provided for symmetric crown types ", paste(symmetric_types, collapse = ", "),
      ". A unique radius will be computed as the mean for tree(s) with id: ",
      paste(trees_inv$id_tree[symmetric_rows][different_radii], collapse = ", "), call. = FALSE
    )
  }
  
  # directional crowns
  directional_rows <- trees_inv$crown_type %in% directional_types
  if (any(directional_rows)) {
    missing_dir_radii <- apply(trees_inv[directional_rows, radii_cols, drop = FALSE], 1, function(x) any(is.na(x)))
    if (any(missing_dir_radii)) stop(
      "All crown radii must be provided for crown types ", paste(directional_types, collapse = ", "),
      ". Missing for tree(s) with id: ", paste(trees_inv$id_tree[directional_rows][missing_dir_radii], collapse = ", "), call. = FALSE
    )
  }
  
  ## ---- conditional hmax_m checks --------------------------------------------
  hmax_types <- c("2E", "8E")
  need_hmax <- any(trees_inv$crown_type %in% hmax_types)
  if ("hmax_m" %in% names(trees_inv)) {
    if (!(is.numeric(trees_inv$hmax_m) || all(is.na(trees_inv$hmax_m)))) stop("`hmax_m` must be numeric or NA.", call. = FALSE)
    if (need_hmax) {
      
      missing_hmax <- is.na(trees_inv$hmax_m) & trees_inv$crown_type %in% hmax_types
      
      if (any(missing_hmax)) stop(
        "`hmax_m` is mandatory for crown_type ", paste(hmax_types, collapse = ", "),
        ". Missing for tree(s) with id: ", paste(trees_inv$id_tree[missing_hmax], collapse = ", "), call. = FALSE
      )
      
      bad_hmax <- trees_inv$hmax_m < trees_inv$hbase_m | trees_inv$hmax_m > trees_inv$h_m
      if (any(bad_hmax, na.rm = TRUE)) {
        stop(
          "`hmax_m` must be between `hbase_m` and `h_m` for tree(s) with id_tree: ",
          paste(trees_inv$id_tree[bad_hmax], collapse = ", "),
          call. = FALSE
        )
      }
      
    } else {
      if (all(is.na(trees_inv$hmax_m))) NULL else if (verbose) warning(
        "`hmax_m` is provided but will be ignored and recomputed.", call. = FALSE
      )
    }
    
  } else if (need_hmax) stop(
    "The column `hmax_m` is required because at least one tree has a crown_type in ", paste(hmax_types, collapse = ", "), call. = FALSE
  )
  
  
  ## ---- height consistency ---------------------------------------------------
  bad_hbase <- trees_inv$hbase_m >= trees_inv$h_m
  if (any(bad_hbase, na.rm = TRUE)) {
    stop(
      "`hbase_m` must be strictly lower than `h_m` for tree(s) with id_tree: ",
      paste(trees_inv$id_tree[bad_hbase], collapse = ", "),
      call. = FALSE
    )
  }
  
  
  ## ---- crown interception properties ---------------------------------------
  has_openness <- "crown_openness" %in% names(trees_inv)
  has_lad <- "crown_lad" %in% names(trees_inv)
  if (!has_openness && !has_lad) stop("Inventory must contain at least one crown interception variable (`crown_openness` or `crown_lad`).", call. = FALSE)
  
  openness_vals <- if (has_openness) trees_inv$crown_openness else rep(NA_real_, nrow(trees_inv))
  lad_vals <- if (has_lad) trees_inv$crown_lad else rep(NA_real_, nrow(trees_inv))
  
  if (has_openness && !(is.numeric(openness_vals) || all(is.na(openness_vals)))) stop("`crown_openness` must be numeric or NA.", call. = FALSE)
  if (has_lad && !(is.numeric(lad_vals) || all(is.na(lad_vals)))) stop("`crown_lad` must be numeric or NA.", call. = FALSE)
  
  missing_interception <- is.na(openness_vals) & is.na(lad_vals)
  if (any(missing_interception)) stop(
    "Each tree must have at least one crown interception property defined (`crown_openness` or `crown_lad`). Missing for tree(s) with id: ", paste(trees_inv$id_tree[missing_interception], collapse = ", "), call. = FALSE
  )
  
  if (has_openness) {
    invalid_openness <- which(trees_inv$crown_openness < 0 | trees_inv$crown_openness > 1)
    if (length(invalid_openness) > 0) stop(
      "`crown_openness` must be between 0 and 1 for tree(s) with id_tree: ",
      paste(trees_inv$id_tree[invalid_openness], collapse = ", "),
      call. = FALSE
    )
  }
  
  if (has_lad) {
    invalid_lad <- which(trees_inv$crown_lad < 0)
    if (length(invalid_lad) > 0) stop(
      "`crown_lad` must be non-negative for tree(s) with id_tree: ",
      paste(trees_inv$id_tree[invalid_lad], collapse = ", "),
      call. = FALSE
    )
  }
  
  if (verbose) {
    if (!has_openness || all(is.na(openness_vals))) warning(
      "`crown_openness` missing or all NA. Porous envelope interception unavailable.", call. = FALSE
    )
    if (!has_lad || all(is.na(lad_vals))) warning(
      "`crown_lad` missing or all NA. Turbid medium interception unavailable.", call. = FALSE
    )
  }
  
  ## ---- success --------------------------------------------------------------
  if (verbose) message("Inventory table successfully validated.")
  invisible(TRUE)
}
