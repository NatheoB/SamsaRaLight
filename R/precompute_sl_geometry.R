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
#'   \item \code{trees}: One row per tree (\code{id_tree}), in the order of
#'     \code{sl_stand$trees}, which is the order of the columns of the
#'     matrices \code{L} and \code{C}.
#'   \item \code{rays}: One row per ray (\code{id_ray}, \code{azimut},
#'     \code{height_angle}, \code{direct}, \code{e_incident} as created for
#'     \link{run_sl_advanced}), with its energy above the canopy on a
#'     horizontal plane (\code{e_horizontal}) and on the slope plane
#'     (\code{e_slope}), in MJ/m2.
#'   \item \code{energies}: Energies above the canopy (direct and diffuse, on
#'     the slope plane and on a horizontal plane), as computed for
#'     \link{run_sl_advanced}.
#'   \item \code{sensors}, \code{cells}: The geometry toward the sensors and
#'     the cells (\code{NULL} if there are no sensors, or if
#'     \code{sensors_only = TRUE} for the cells). Energies are on a horizontal
#'     plane for sensors and on the slope plane for cells, as in
#'     \link{run_sl_advanced}. Each one is a list with:
#'     \itemize{
#'       \item \code{targets}: One row per target (\code{id_sensor} or
#'         \code{id_cell}), with the energy of the rays reaching the target
#'         without crossing any crown or trunk (\code{e_unobs_direct},
#'         \code{e_unobs_diffuse}, in MJ/m2).
#'       \item \code{pairs}: One row per ray x target pair intercepted at
#'         least once by a crown, and not by a trunk, with the index of the
#'         target (\code{target}, row in \code{targets}) and of the ray
#'         (\code{ray}, row in \code{rays}), and the number of crowns
#'         crossed (\code{n_trees}).
#'       \item \code{L}: Sparse matrix (pairs x trees) of the path length of
#'         the ray of each pair across each crown (in m, summed over crown
#'         parts and torus copies of the tree).
#'       \item \code{C}: Sparse pattern matrix (pairs x trees) of the crowns
#'         crossed by the ray of each pair (a crown is crossed even if the
#'         path length is null, e.g. a ray tangent to the crown, which matters
#'         for porous envelopes).
#'       \item \code{W}: Sparse matrix (targets x pairs) of the energy above
#'         the canopy of the ray of each pair (in MJ/m2), i.e. the energy that
#'         would reach the target if the crossed crowns were transparent.
#'       \item \code{trunks}: The same \code{pairs}, \code{L}, \code{C} and
#'         \code{W} for the ray x target pairs intercepted by a trunk (with
#'         \code{n_trees = 0} if no crown is crossed). They bring no energy if
#'         \code{trunk_interception = TRUE} in
#'         \link{estimate_pacl_from_geometry}, and are kept apart so that
#'         they are not computed in this case.
#'     }
#'   \item \code{params}: The settings used to compute the geometry, and the
#'     version of the format of the object (\code{format_version}).
#' }
#' The object should not be modified: filter the output of
#' \link{estimate_pacl_from_geometry} instead.
#'
#' @seealso \link{estimate_pacl_from_geometry}, \link{run_sl_advanced}
#'
#' @importFrom Matrix sparseMatrix
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

  # Geometry toward a set of targets as sparse matrices (NULL if not computed),
  # with the pairs intercepted by a trunk kept apart (in trunks)
  # crossings: one row per crown crossed by the ray of a pair (pair, tree, length), as exported by the C++
  n_trees <- nrow(sl_stand$trees)
  as_matrices <- function(geom, rays_e) {

    if (is.null(geom)) return(NULL)

    n_targets <- nrow(geom$targets)
    pairs <- geom$pairs
    crossings <- geom$crossings

    # Pairs selected by keep, with their crossings (pairs renumbered within the selection)
    select_pairs <- function(keep) {

      n_pairs <- sum(keep)
      index <- cumsum(keep)
      keep_crossings <- keep[crossings$pair]

      list(
        pairs = data.frame(target = pairs$target[keep], ray = pairs$ray[keep], n_trees = pairs$n_trees[keep]),
        L = Matrix::sparseMatrix(i = index[crossings$pair[keep_crossings]], j = crossings$tree[keep_crossings],
                                 x = crossings$length[keep_crossings], dims = c(n_pairs, n_trees)),
        C = Matrix::sparseMatrix(i = index[crossings$pair[keep_crossings]], j = crossings$tree[keep_crossings],
                                 dims = c(n_pairs, n_trees)),
        W = Matrix::sparseMatrix(i = pairs$target[keep], j = seq_len(n_pairs), x = rays_e[pairs$ray[keep]],
                                 dims = c(n_targets, n_pairs))
      )
    }

    c(list(targets = geom$targets), select_pairs(!pairs$blocked), list(trunks = select_pairs(pairs$blocked)))
  }

  # Rays as created for run_sl_advanced(), with their energy above the canopy on both planes (computed in C++)
  rays <- monthly_rays$rays[, c("id_ray", "azimut", "height_angle", "direct", "e_incident")]
  rays$e_horizontal <- out$rays$e_horizontal
  rays$e_slope <- out$rays$e_slope

  sl_geometry <- list(
    trees = data.frame(id_tree = sl_stand$trees$id_tree),
    rays = rays,
    energies = monthly_rays$energies,
    sensors = as_matrices(out$sensors, rays$e_horizontal),
    cells = as_matrices(out$cells, rays$e_slope),
    params = list(
      "format_version" = 2L,
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
#' Checks the structure of an object created by \link{precompute_sl_geometry}
#' (format version 2): its components, the links between \code{targets},
#' \code{pairs}, \code{rays} and \code{trees}, and the dimensions and
#' content of the matrices \code{L}, \code{C} and \code{W}, for the pairs
#' not intercepted and intercepted by a trunk.
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

  if (!identical(x$params$format_version, 2L)) {
    stop("`sl_geometry` was created by an older version of SamsaRaLight: ",
         "compute it again with precompute_sl_geometry().", call. = FALSE)
  }

  required <- c("trees", "rays", "energies", "sensors", "cells", "params")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Missing element(s) in `sl_geometry`: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (is.null(x$sensors) && is.null(x$cells)) {
    stop("`sl_geometry` must contain the geometry toward sensors or cells.", call. = FALSE)
  }

  n_rays <- nrow(x$rays)
  n_trees <- nrow(x$trees)

  # ---- Links between the tables and the matrices, for sensors and cells ----
  for (type in c("sensors", "cells")) {

    geom <- x[[type]]
    if (is.null(geom)) next

    n_targets <- nrow(geom$targets)

    # Pairs not intercepted (geom) and intercepted (geom$trunks) by a trunk
    for (block in c("pairs", "trunks")) {

      b <- if (block == "pairs") geom else geom$trunks
      n_pairs <- nrow(b$pairs)
      where <- if (block == "pairs") type else paste0(type, "$trunks")
      err <- function(msg) stop(sprintf("%s in `sl_geometry$%s`.", msg, where), call. = FALSE)

      # Pairs
      if (any(b$pairs$target < 1 | b$pairs$target > n_targets)) err("Indices of targets out of range in `pairs`")
      if (any(b$pairs$ray < 1 | b$pairs$ray > n_rays)) err("Indices of rays out of range in `pairs`")

      # Matrices: classes and dimensions
      if (!inherits(b$L, "dgCMatrix") || !identical(dim(b$L), c(n_pairs, n_trees))) {
        err("`L` must be a sparse matrix (dgCMatrix) of dimensions pairs x trees")
      }
      if (!inherits(b$C, "ngCMatrix") || !identical(dim(b$C), c(n_pairs, n_trees))) {
        err("`C` must be a sparse pattern matrix (ngCMatrix) of dimensions pairs x trees")
      }
      if (!inherits(b$W, "dgCMatrix") || !identical(dim(b$W), c(n_targets, n_pairs))) {
        err("`W` must be a sparse matrix (dgCMatrix) of dimensions targets x pairs")
      }

      # Content: crowns crossed per pair, non-negative path lengths, one target per pair (column of W)
      if (!identical(as.integer(Matrix::rowSums(b$C)), b$pairs$n_trees)) {
        err("Number of crowns crossed per pair inconsistent between `pairs` and `C`")
      }
      if (any(b$L@x < 0)) err("Negative path lengths in `L`")
      if (!identical(diff(b$W@p), rep(1L, n_pairs)) || !identical(b$W@i + 1L, b$pairs$target)) {
        err("Each column of `W` must have a single value, in the row of the target of the pair")
      }
    }

    # Pairs not intercepted by a trunk cross at least one crown (otherwise they would be unobstructed rays)
    if (any(geom$pairs$n_trees == 0)) {
      stop(sprintf("Pairs crossing no crown and no trunk in `sl_geometry$%s$pairs`.", type), call. = FALSE)
    }
  }

  invisible(TRUE)
}
