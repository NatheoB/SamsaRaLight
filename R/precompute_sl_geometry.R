#' Precompute the geometry of light interception in a stand
#'
#' Runs, once for a given stand, the part of the SamsaRaLight ray tracing that
#' does not depend on the crown leaf area density (LAD) or openness: for each
#' ray toward each target (sensor or cell), the crowns crossed by the ray, the
#' path length across each crown, and whether a trunk intercepts the ray.
#' The light reaching the targets can then be estimated very quickly for any
#' crown LAD or openness with \link{estimate_pacl_from_geometry}, e.g. when
#' calibrating crown LAD from light measurements, where the same stand is run
#' many times with different LAD values.
#'
#' @inheritParams run_sl_advanced
#'
#' @details
#' The energy reaching a target does not depend on the order in which a ray
#' crosses the crowns: it is the energy of the ray above the canopy times the
#' product of the transmittances of the crossed crowns (Beer-Lambert law for a
#' turbid medium, crown openness for a porous envelope), and is null if a
#' trunk intercepts the ray. The light on sensors and cells estimated by
#' \link{estimate_pacl_from_geometry} is thus the same as with
#' \link{run_sl_advanced} (up to rounding errors) for the same stand,
#' radiations and ray settings. The light intercepted by trees (\code{e},
#' \code{epot}, \code{lci}, ...) depends on the order of interceptions and
#' cannot be estimated from the geometry.
#'
#' The geometry depends on the stand (tree positions, dimensions and crown
#' shapes, sensors, cells), the radiations and the ray settings, but not on
#' the crown LAD and openness, the extinction coefficient, the clumping factor
#' or the interception by trunks, which are chosen in
#' \link{estimate_pacl_from_geometry}. The geometry must be computed again if
#' the stand, the radiations or the ray settings change. Trunks are always
#' included in the geometry.
#'
#' The size of the geometry increases with the number of targets, rays and
#' crowns crossed by each ray: use \code{sensors_only = TRUE} if the light on
#' the cells is not needed.
#'
#' @return An object of class \code{"sl_geometry"} (a list), to be used in
#'   \link{estimate_pacl_from_geometry}, and that can be saved with
#'   \code{saveRDS()}. It contains:
#' \itemize{
#'   \item \code{sensors}, \code{cells}: The geometry toward the sensors and
#'     the cells (\code{NULL} if there are no sensors, or if
#'     \code{sensors_only = TRUE} for the cells). Each one is a list with:
#'     \itemize{
#'       \item \code{targets}: One row per target (\code{id_sensor} or
#'         \code{id_cell}), with the energy of the rays reaching the target
#'         without crossing any crown or trunk (\code{e_unobs_direct},
#'         \code{e_unobs_diffuse}, in MJ/m2 on a horizontal plane for sensors
#'         and on the slope plane for cells, as in \link{run_sl_advanced}).
#'       \item \code{pairs}: One row per ray x target pair intercepted at
#'         least once by a crown or a trunk, with the index of the target
#'         (\code{target}, row in \code{targets}) and of the ray (\code{ray},
#'         row in \code{rays}), whether a trunk intercepts the ray
#'         (\code{blocked}), and the number of crowns crossed
#'         (\code{n_trees}).
#'       \item \code{crossings}: One row per crown crossed by the ray of a
#'         pair, with the index of the pair (\code{pair}, row in
#'         \code{pairs}), of the tree (\code{tree}, row in
#'         \code{sl_stand$trees}) and the path length across the crown
#'         (\code{length}, in m, summed over crown parts).
#'     }
#'   \item \code{rays}: One row per ray, with its energy above the canopy
#'     (in MJ/m2, on a horizontal plane \code{e_horizontal} and on the slope
#'     plane \code{e_slope}) and its type (\code{direct}).
#'   \item \code{id_tree}: Identifiers of the trees, in the order of
#'     \code{sl_stand$trees} (order of the crown LAD or openness to give in
#'     \link{estimate_pacl_from_geometry}).
#'   \item \code{energies}: Energies above the canopy, as computed for
#'     \link{run_sl_advanced}.
#'   \item \code{params}: The settings used to compute the geometry.
#' }
#' The object should not be modified: filter the output of
#' \link{estimate_pacl_from_geometry} instead.
#'
#' @seealso \link{estimate_pacl_from_geometry}, \link{run_sl_advanced}
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
#' str(geometry, max.level = 2)
#'
#' @export
#'
precompute_sl_geometry <- function(
    sl_stand,
    monthly_radiations,
    sensors_only = FALSE,
    use_torus = TRUE,
    height_anglemin = 10,
    direct_startoffset = 0,
    direct_anglestep = 5,
    diffuse_anglestep = 15,
    soc = TRUE,
    start_day = 1,
    end_day = 365,
    parallel_mode = FALSE,
    n_threads = NULL,
    verbose = TRUE
) {

  # ---- ARGUMENTS CHECKS ----

  check_sl_run_args(
    sl_stand,
    monthly_radiations,
    logical_args = list(
      sensors_only = sensors_only, use_torus = use_torus, soc = soc,
      parallel_mode = parallel_mode, verbose = verbose
    ),
    numeric_args = list(
      height_anglemin = height_anglemin, direct_startoffset = direct_startoffset,
      direct_anglestep = direct_anglestep, diffuse_anglestep = diffuse_anglestep
    ),
    start_day = start_day,
    end_day = end_day,
    n_threads = n_threads
  )

  # Targets
  if (sensors_only && NROW(sl_stand$sensors) == 0) {
    stop("`sensors_only = TRUE` requires sensors in `sl_stand`.", call. = FALSE)
  }


  # CREATE RAYS ----
  monthly_rays <- create_sl_rays(
    monthly_rad = monthly_radiations,
    latitude = sl_stand$geometry$latitude,
    start_day = start_day,
    end_day = end_day,
    soc = soc,
    slope = sl_stand$geometry$slope,
    north_to_x_cw = sl_stand$geometry$north2x,
    aspect = sl_stand$geometry$aspect,
    height_anglemin = height_anglemin,
    direct_startoffset = direct_startoffset,
    direct_anglestep = direct_anglestep,
    diffuse_anglestep = diffuse_anglestep
  )


  # RUN C++ GEOMETRY ----

  ## Set parallel (BLAS single-threaded, OpenMP threads) ----
  set_sl_parallel(parallel_mode, n_threads, verbose)

  ## Compute the geometry ----
  out <- precompute_sl_geometry_rcpp(
    sl_stand$trees,
    sl_stand$sensors,
    sensors_only,
    sl_stand$cells,
    monthly_rays$rays,
    monthly_rays$energies[["slope_direct"]],
    monthly_rays$energies[["slope_diffuse"]],
    monthly_rays$energies[["horizontal_direct"]],
    monthly_rays$energies[["horizontal_diffuse"]],
    sl_stand$geometry$slope,
    sl_stand$geometry$north2x,
    sl_stand$geometry$aspect,
    sl_stand$geometry$cell_size,
    sl_stand$geometry$n_cells_x,
    sl_stand$geometry$n_cells_y,
    use_torus
  )


  # PREPARE OUTPUT ----

  # Interceptions with an odd number of interception points that are not a known particular case
  # (ignored as in Capsis, but reported to the user)
  warn_sl_unexpected_interceptions(out$n_unexpected_interceptions)

  sl_geometry <- list(
    sensors = out$sensors,
    cells = out$cells,
    rays = out$rays,
    id_tree = sl_stand$trees$id_tree,
    energies = monthly_rays$energies,
    params = list(
      "sensors_only" = sensors_only,
      "start_day" = start_day,
      "end_day" = end_day,
      "soc" = soc,
      "use_torus" = use_torus,
      "height_anglemin" = height_anglemin,
      "direct_startoffset" = direct_startoffset,
      "direct_anglestep" = direct_anglestep,
      "diffuse_anglestep" = diffuse_anglestep
    )
  )

  class(sl_geometry) <- c("sl_geometry", "list")

  validate_sl_geometry(sl_geometry)
  if (verbose) message("SamsaRaLight geometry was computed successfully.")

  return(sl_geometry)
}


#' Validate an sl_geometry object
#'
#' Checks the structure of an object created by \link{precompute_sl_geometry},
#' and the links between its tables (indices of targets, rays, pairs and
#' trees, and number of crowns crossed per pair).
#'
#' @param x An object of class \code{"sl_geometry"}.
#'
#' @return Invisibly \code{TRUE}, or an error.
#'
#' @keywords internal
#'
validate_sl_geometry <- function(x) {

  # ---- Class and components ----
  if (!inherits(x, "sl_geometry")) {
    stop("Object must inherit from class 'sl_geometry'.", call. = FALSE)
  }

  required <- c("sensors", "cells", "rays", "id_tree", "energies", "params")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Missing element(s) in `sl_geometry`: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (is.null(x$sensors) && is.null(x$cells)) {
    stop("`sl_geometry` must contain the geometry toward sensors or cells.", call. = FALSE)
  }

  n_rays <- nrow(x$rays)
  n_trees <- length(x$id_tree)

  # ---- Links between the tables, for sensors and cells ----
  for (type in c("sensors", "cells")) {

    geom <- x[[type]]
    if (is.null(geom)) next

    n_targets <- nrow(geom$targets)
    n_pairs <- nrow(geom$pairs)

    if (any(geom$pairs$target < 1 | geom$pairs$target > n_targets)) {
      stop(sprintf("Indices of targets out of range in `sl_geometry$%s$pairs`.", type), call. = FALSE)
    }
    if (any(geom$pairs$ray < 1 | geom$pairs$ray > n_rays)) {
      stop(sprintf("Indices of rays out of range in `sl_geometry$%s$pairs`.", type), call. = FALSE)
    }
    if (any(geom$crossings$pair < 1 | geom$crossings$pair > n_pairs)) {
      stop(sprintf("Indices of pairs out of range in `sl_geometry$%s$crossings`.", type), call. = FALSE)
    }
    if (any(geom$crossings$tree < 1 | geom$crossings$tree > n_trees)) {
      stop(sprintf("Indices of trees out of range in `sl_geometry$%s$crossings`.", type), call. = FALSE)
    }
    if (!identical(tabulate(geom$crossings$pair, nbins = n_pairs), geom$pairs$n_trees)) {
      stop(sprintf("Number of crowns crossed per pair inconsistent between `pairs` and `crossings` in `sl_geometry$%s`.", type),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}
