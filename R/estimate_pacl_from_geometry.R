#' Estimate the light on sensors and cells from a precomputed geometry
#'
#' Estimates the light reaching the sensors and cells of a stand for given
#' crown leaf area densities (turbid medium) or crown openness (porous
#' envelope), from the geometry computed by \link{precompute_sl_geometry}.
#' Results are the same as the \code{sensors} and \code{cells} outputs of
#' \link{run_sl_advanced} with the same settings (up to rounding errors), but
#' much faster, since the rays are not traced again. The light intercepted by
#' trees cannot be estimated from the geometry.
#'
#' @inheritParams run_sl_advanced
#' @param sl_geometry An object of class \code{"sl_geometry"} created by
#'   \link{precompute_sl_geometry}.
#' @param trees A data frame with one row per tree of the stand used to
#'   compute the geometry, with the column \code{id_tree}, and the column
#'   \code{crown_lad} (if \code{turbid_medium = TRUE}) or
#'   \code{crown_openness} (if \code{turbid_medium = FALSE}), as in
#'   \code{sl_stand$trees}. Other columns are ignored, and the order of the
#'   rows does not matter (trees are matched by \code{id_tree}).
#'
#' @details
#' For each ray toward each target, the energy reaching the target is the
#' energy of the ray above the canopy times the transmittance of the crossed
#' crowns: \eqn{\exp(-k \Omega \sum LAD \times length)} for turbid media, or
#' the product of the crown openness for porous envelopes. It is null if a
#' trunk intercepts the ray and \code{trunk_interception = TRUE}.
#'
#' @return A list with two data frames, \code{sensors} and \code{cells}, with
#'   the same columns as in the output of \link{run_sl_advanced}
#'   (\code{light$sensors} and \code{light$cells}): \code{id_sensor} or
#'   \code{id_cell}, \code{e}, \code{pacl} and \code{punobs}, and their direct
#'   and diffuse components if \code{detailed_output = TRUE}. As in
#'   \link{run_sl_advanced}, energies are on a horizontal plane for sensors and
#'   on the slope plane for cells. A data frame is empty if the geometry toward
#'   these targets was not computed.
#'
#' @seealso \link{precompute_sl_geometry}, \link{run_sl_advanced}
#'
#' @examples
#' data_prenovel <- SamsaRaLight::data_prenovel
#'
#' stand <- create_sl_stand(
#'   trees_inv = data_prenovel$trees,
#'   cell_size = 10,
#'   latitude = data_prenovel$info$latitude,
#'   slope = data_prenovel$info$slope,
#'   aspect = data_prenovel$info$aspect,
#'   north2x = data_prenovel$info$north2x
#' )
#'
#' geometry <- precompute_sl_geometry(
#'   sl_stand = stand,
#'   monthly_radiations = data_prenovel$radiations
#' )
#'
#' # Light on the cells for two crown LAD values
#' trees <- stand$trees
#'
#' trees$crown_lad <- 0.5
#' head(estimate_pacl_from_geometry(geometry, trees)$cells)
#'
#' trees$crown_lad <- 1
#' head(estimate_pacl_from_geometry(geometry, trees)$cells)
#'
#' @export
#'
estimate_pacl_from_geometry <- function(
    sl_geometry,
    trees,
    turbid_medium = TRUE,
    extinction_coef = 0.5,
    clumping_factor = 1,
    trunk_interception = TRUE,
    detailed_output = FALSE
) {

  # ---- ARGUMENTS CHECKS ----
  # Light checks, as the function is meant to be called many times
  # (the structure of the geometry is checked when it is created)

  if (!inherits(sl_geometry, "sl_geometry")) {
    stop("`sl_geometry` must be an object created by precompute_sl_geometry().", call. = FALSE)
  }

  for (arg_name in c("turbid_medium", "trunk_interception", "detailed_output")) {
    arg_val <- get(arg_name)
    if (!is.logical(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single TRUE or FALSE.", arg_name), call. = FALSE)
    }
  }

  for (arg_name in c("extinction_coef", "clumping_factor")) {
    arg_val <- get(arg_name)
    if (!is.numeric(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single numeric value.", arg_name), call. = FALSE)
    }
  }

  # Crown LAD or openness (same checks as in run_sl_advanced())
  validate_interception_model(trees, turbid_medium)

  # Trees matched by id_tree with the trees of the geometry
  idx <- match(sl_geometry$id_tree, trees$id_tree)
  if (anyNA(idx) || anyDuplicated(trees$id_tree) || nrow(trees) != length(sl_geometry$id_tree)) {
    stop("`trees$id_tree` must contain exactly the trees of the stand used to compute the geometry.",
         call. = FALSE)
  }

  # Attenuation coefficient per unit path length in each crown (turbid medium) or crown openness (porous envelope)
  if (turbid_medium) {
    tree_coef <- extinction_coef * clumping_factor * trees$crown_lad[idx]
  } else {
    tree_coef <- trees$crown_openness[idx]
  }


  # ---- ENERGIES ON THE TARGETS ----
  # As in the C++ export of run_sl_advanced(): energy above canopy minus intercepted energy (MJ/m2),
  # on a horizontal plane for sensors, on the slope plane for cells
  estimate_targets <- function(geom, name_id, rays_e, e_above_direct, e_above_diffuse) {

    # Targets not computed in the geometry: no target
    if (is.null(geom)) {
      ids <- integer(0)
      e_direct <- e_diffuse <- e_unobs_direct <- e_unobs_diffuse <- numeric(0)
    } else {

      energies <- estimate_pacl_from_geometry_rcpp(
        geom$targets$e_unobs_direct, geom$targets$e_unobs_diffuse,
        geom$pairs$target, geom$pairs$ray, geom$pairs$blocked, geom$pairs$n_trees,
        geom$crossings$pair, geom$crossings$tree, geom$crossings$length,
        rays_e, sl_geometry$rays$direct,
        as.numeric(tree_coef), turbid_medium, trunk_interception
      )

      ids <- geom$targets[[name_id]]
      e_direct <- e_above_direct - energies[, "e_intercepted_direct"]
      e_diffuse <- e_above_diffuse - energies[, "e_intercepted_diffuse"]
      e_unobs_direct <- energies[, "e_unobs_direct"]
      e_unobs_diffuse <- energies[, "e_unobs_diffuse"]
    }

    e <- e_direct + e_diffuse
    pacl <- e / (e_above_direct + e_above_diffuse)
    punobs <- (e_unobs_direct + e_unobs_diffuse) / e

    # Columns of the output of run_sl_advanced(), built as a list of vectors then converted with list2DF()
    # (much faster than data.frame() when the function is called many times, e.g. in a calibration)
    if (detailed_output) {
      out <- list(ids, e, e_direct, e_diffuse,
                  pacl, e_direct / e_above_direct, e_diffuse / e_above_diffuse,
                  punobs, e_unobs_direct / e_direct, e_unobs_diffuse / e_diffuse)
      names(out) <- c(name_id, "e", "e_direct", "e_diffuse",
                      "pacl", "pacl_direct", "pacl_diffuse",
                      "punobs", "punobs_direct", "punobs_diffuse")
    } else {
      out <- list(ids, e, pacl, punobs)
      names(out) <- c(name_id, "e", "pacl", "punobs")
    }

    list2DF(out)
  }

  list(
    sensors = estimate_targets(sl_geometry$sensors, "id_sensor", sl_geometry$rays$e_horizontal,
                               sl_geometry$energies[["horizontal_direct"]], sl_geometry$energies[["horizontal_diffuse"]]),
    cells = estimate_targets(sl_geometry$cells, "id_cell", sl_geometry$rays$e_slope,
                             sl_geometry$energies[["slope_direct"]], sl_geometry$energies[["slope_diffuse"]])
  )
}
