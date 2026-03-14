#' Compute advanced SamsaRaLight radiative balance
#'
#' This function runs the full light interception and radiative balance
#' simulation for a virtual forest stand with advanced parameters. It allows
#' customization of ray discretization, sky type and trunk interception.  
#'
#' For typical use, see the simpler \link{run_sl} wrapper that sets standard
#' discretization parameters for most users.
#'
#' @param sl_stand An object of class \code{"sl_stand"} representing the virtual stand.
#'   Each row is a tree with required and optional columns describing crown geometry,
#'   height, crown radius, crown openness, LAD, etc. See \link{validate_sl_stand}.
#' @param monthly_radiations data.frame of monthly horizontal radiation (Hrad) and
#'   diffuse to global ratio (DGratio), computed with \link{get_monthly_radiations}.
#' @param sensors_only logical, if TRUE, compute interception only for sensors
#' @param use_torus logical, if TRUE, use torus system for borders
#' @param turbid_medium logical, if TRUE, crowns are considered turbid medium (using column `crown_lad`), else porous envelope (using column `crown_openess`)
#' @param extinction_coef Numeric scalar. Leaf extinction coefficient controlling
#'   the probability that a ray is intercepted by foliage. It represents the
#'   effective light attenuation per unit leaf area and is linked to average
#'   leaf orientation. Higher values increase interception (default = 0.5).
#' @param clumping_factor Numeric scalar controlling the aggregation of leaves
#'   within the crown volume. A value of 1 corresponds to a homogeneous (random)
#'   foliage distribution; values < 1 indicate clumped foliage, and values > 1
#'   indicate more regular spacing. This modifies effective light interception
#'   in the turbid medium model (default = 1).
#' @param trunk_interception logical, if TRUE, account for trunk interception
#' @param height_anglemin numeric, minimum altitude angle for rays (degrees)
#' @param direct_startoffset numeric, starting angle of first direct ray (degrees)
#' @param direct_anglestep numeric, hour angle step between direct rays (degrees)
#' @param diffuse_anglestep numeric, hour angle step between diffuse rays (degrees)
#' @param soc logical, if TRUE, use Standard Overcast Sky; if FALSE, Uniform Overcast Sky
#' @param start_day integer, first day of the vegetative period (1–365)
#' @param end_day integer, last day of the vegetative period (1–365)
#' @param detailed_output logical, if TRUE, include detailed rays, energies, and interception matrices
#' @param parallel_mode logical. If TRUE, ray–target computations are parallelised
#'   using OpenMP. If FALSE, the model runs in single-thread mode.
#' @param n_threads integer or NULL. Number of CPU threads to use when
#'   \code{parallel_mode = TRUE}. If NULL (default), OpenMP automatically selects
#'   the number of available cores. If provided, must be a positive integer.
#' @param verbose Logical; if \code{TRUE}, informative messages are printed.
#'
#' @return An object of class \code{"sl_output"} (list) containing:
#' \itemize{
#'   \item \code{light}: list with simulation outputs for trees, cells, and sensors
#'   \item \code{info}: list with run metadata (latitude, days, sky type, etc.)
#'   \item \code{monthly_rays} (if detailed_output = TRUE): ray discretization per month
#'   \item \code{interceptions} (if detailed_output = TRUE): tree/cell interception matrices
#' }
#'
#' @details
#' This advanced function exposes all ray tracing parameters and is intended
#' for users who need full control over ray discretization and modeling options.
#' For most users, see \link{run_sl} which wraps this function with default
#' parameters suitable for standard runs.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr select %>%
#' 
#' @useDynLib SamsaRaLight, .registration = TRUE
#' 
#' @export
#'
run_sl_advanced <- function(
    sl_stand,
    monthly_radiations,
    sensors_only = FALSE,
    use_torus = TRUE,
    turbid_medium = TRUE,
    extinction_coef = 0.5,
    clumping_factor = 1,
    trunk_interception = TRUE,
    height_anglemin = 10,
    direct_startoffset = 0,
    direct_anglestep = 5,
    diffuse_anglestep = 15,
    soc = TRUE,
    start_day = 1,
    end_day = 365,
    detailed_output = FALSE,
    parallel_mode = FALSE,
    n_threads = NULL,
    verbose = TRUE
) {
  
  # ---- ARGUMENTS CHECKS ----
  
  # sl_stand
  validate_sl_stand(sl_stand)
  
  # monthly_radiations
  check_monthly_radiations(monthly_radiations, verbose = FALSE)
  
  # Logical arguments
  for (arg_name in c("sensors_only", "use_torus", "turbid_medium", "trunk_interception",
                     "soc", "detailed_output", "parallel_mode", "verbose")) {
    arg_val <- get(arg_name)
    if (!is.logical(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single TRUE or FALSE.", arg_name), call. = FALSE)
    }
  }
  
  # Numeric scalar arguments
  for (arg_name in c("extinction_coef", "clumping_factor",
                     "height_anglemin", "direct_startoffset",
                     "direct_anglestep", "diffuse_anglestep")) {
    arg_val <- get(arg_name)
    if (!is.numeric(arg_val) || length(arg_val) != 1 || is.na(arg_val)) {
      stop(sprintf("`%s` must be a single numeric value.", arg_name), call. = FALSE)
    }
  }
  
  # Days
  if (!is.numeric(start_day) || length(start_day) != 1 || start_day < 1 || start_day > 365) {
    stop("`start_day` must be a single number between 1 and 365.", call. = FALSE)
  }
  if (!is.numeric(end_day) || length(end_day) != 1 || end_day < 1 || end_day > 365) {
    stop("`end_day` must be a single number between 1 and 365.", call. = FALSE)
  }
  if (end_day < start_day) {
    stop("`end_day` cannot be smaller than `start_day`.", call. = FALSE)
  }
  
  # Parallel threads
  if (!is.null(n_threads)) {
    if (!is.numeric(n_threads) || length(n_threads) != 1 || is.na(n_threads)) {
      stop("`n_threads` must be a single integer or NULL.", call. = FALSE)
    }
    if (n_threads < 1 || n_threads %% 1 != 0) {
      stop("`n_threads` must be a positive integer.", call. = FALSE)
    }
  }
  
  # Interception model checks
  validate_interception_model(sl_stand$trees, turbid_medium)
  
  
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
  
  # RUN C++ SIMULATION ----

  ## Prevent BLAS and OpenMP fighting ----
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)

  ## Set parallel ----
  if (is.null(n_threads)) n_threads <- -1L
  sl_set_openmp(parallel_mode = parallel_mode, 
                num_threads = as.integer(n_threads),
                verbose = verbose)
  
  if (verbose) sl_print_openmp_status()
  
  
  ## Run the model ----
  out <- sl_run_rcpp(
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
    use_torus, 
    turbid_medium, 
    extinction_coef,
    clumping_factor,
    trunk_interception
  )
  
  # PREPARE OUTPUT ----
  interceptions <- out$interceptions
  out$interceptions <- NULL
  
  if (!detailed_output) {
    out$sensors <- out$sensors %>% dplyr::select(id_sensor, e, pacl, punobs)
    out$cells   <- out$cells %>% dplyr::select(id_cell, e, pacl, punobs)
    out$trees   <- out$trees %>% dplyr::select(id_tree, epot, e, lci, eunobs)
  }
  
  out_sl <- list(
    output = list(
      "light" = out
    ),
    input = list(
      "sl_stand" = sl_stand,
      "monthly_radiations" = monthly_radiations,
      "params" = list(
        "detailed_output" = detailed_output,
        "start_day" = start_day,
        "end_day" = end_day,
        "soc" = soc,
        "use_torus" = use_torus,
        "turbid_medium" = turbid_medium,
        "extinction_coef" = extinction_coef,
        "clumping_factor" = clumping_factor,
        "trunk_interception" = trunk_interception,
        "height_anglemin" = height_anglemin,
        "direct_startoffset" = direct_startoffset,
        "direct_anglestep" = direct_anglestep,
        "diffuse_anglestep" = diffuse_anglestep
      )
    )
  )
  
  if (detailed_output) {
    out_sl$output$monthly_rays <- monthly_rays
    out_sl$output$interceptions <- interceptions
  }
  
  class(out_sl) <- c("sl_output", "list")
  
  validate_sl_output(out_sl)
  if (verbose) message("SamsaRaLight simulation was run successfully.")
  
  return(out_sl)
}


#' Validate a SamsaRaLight simulation output object
#'
#' Performs structural and internal consistency checks on an object
#' returned by \code{run_sl_advanced()} or \code{run_sl()}.
#'
#' This function is called internally at the end of a simulation to ensure
#' that the returned object is valid and complete before being provided
#' to the user.
#'
#' @param x Object expected to inherit from class \code{"sl_output"}.
#'
#' @details
#' The function checks that:
#' \itemize{
#'   \item The object inherits from class \code{"sl_output"}.
#'   \item Top-level components \code{$output} and \code{$input} exist.
#'   \item \code{$output$light} contains \code{trees}, \code{cells}, and \code{sensors}.
#'   \item These elements are data.frames.
#'   \item Required identifier and energy columns are present.
#'   \item Energy-related columns are numeric.
#' }
#'
#' If detailed outputs were requested, the presence and structure of
#' \code{$output$monthly_rays} and \code{$output$interceptions}
#' are also verified when available.
#'
#' @return Invisibly returns \code{TRUE} if validation passes.
#'   Stops with an informative error message otherwise.
#'
#' @keywords internal
validate_sl_output <- function(x) {
  
  # ---- Class ----
  if (!inherits(x, "sl_output")) {
    stop("Object must inherit from class 'sl_output'.", call. = FALSE)
  }
  
  # ---- Top-level structure ----
  if (!is.list(x$output)) {
    stop("`x$output` must be a list.", call. = FALSE)
  }
  
  if (!is.list(x$input)) {
    stop("`x$input` must be a list.", call. = FALSE)
  }
  
  if (!"light" %in% names(x$output)) {
    stop("`x$output` must contain element `light`.", call. = FALSE)
  }
  
  light <- x$output$light
  
  # ---- Required light components ----
  required_light <- c("trees", "cells", "sensors")
  missing_light <- setdiff(required_light, names(light))
  if (length(missing_light) > 0) {
    stop("Missing element(s) in `output$light`: ",
         paste(missing_light, collapse = ", "),
         call. = FALSE)
  }
  
  # ---- Check data.frame structure ----
  for (comp in required_light) {
    if (!is.data.frame(light[[comp]])) {
      stop(sprintf("`output$light$%s` must be a data.frame.", comp),
           call. = FALSE)
    }
  }
  
  # ---- Minimal required columns ----
  
  # Trees
  trees_req <- c("id_tree", "e")
  missing_trees <- setdiff(trees_req, names(light$trees))
  if (length(missing_trees) > 0) {
    stop("`trees` output is missing column(s): ",
         paste(missing_trees, collapse = ", "),
         call. = FALSE)
  }
  
  # Cells
  cells_req <- c("id_cell", "e")
  missing_cells <- setdiff(cells_req, names(light$cells))
  if (length(missing_cells) > 0) {
    stop("`cells` output is missing column(s): ",
         paste(missing_cells, collapse = ", "),
         call. = FALSE)
  }
  
  # Sensors
  sensors_req <- c("id_sensor", "e")
  missing_sensors <- setdiff(sensors_req, names(light$sensors))
  if (length(missing_sensors) > 0) {
    stop("`sensors` output is missing column(s): ",
         paste(missing_sensors, collapse = ", "),
         call. = FALSE)
  }
  
  # ---- Energy columns must be numeric ----
  numeric_check <- function(df, cols, name) {
    for (col in intersect(cols, names(df))) {
      if (!is.numeric(df[[col]])) {
        stop(sprintf("Column `%s` in `%s` must be numeric.", col, name),
             call. = FALSE)
      }
    }
  }
  
  numeric_check(light$trees,   c("e", "epot", "lci", "eunobs"), "trees")
  numeric_check(light$cells,   c("e", "pacl", "punobs"), "cells")
  numeric_check(light$sensors, c("e", "pacl", "punobs"), "sensors")
  
  # ---- Optional detailed outputs ----
  if ("monthly_rays" %in% names(x$output)) {
    if (!is.list(x$output$monthly_rays)) {
      stop("`output$monthly_rays` must be a list.", call. = FALSE)
    }
  }
  
  if ("interceptions" %in% names(x$output)) {
    if (!is.list(x$output$interceptions)) {
      stop("`output$interceptions` must be a list.", call. = FALSE)
    }
  }
  
  invisible(TRUE)
}



#' Validate tree interception model parameters against model type
#'
#' @keywords internal
validate_interception_model <- function(trees, turbid_medium) {
  
  if (!is.data.frame(trees)) {
    stop("`trees` must be a data.frame.", call. = FALSE)
  }
  
  if (!"id_tree" %in% names(trees)) {
    stop("`trees` must contain column `id_tree`.", call. = FALSE)
  }
  
  if (turbid_medium) {
    
    if (!"crown_lad" %in% names(trees)) {
      stop(
        "turbid_medium = TRUE requires column `crown_lad` in tree inventory.",
        call. = FALSE
      )
    }
    
    if (any(is.na(trees$crown_lad))) {
      stop(
        "turbid_medium = TRUE requires to define `crown_lad` for all trees.", 
        call. = FALSE)
    }
    
    if (!is.numeric(trees$crown_lad)) {
      stop("turbid_medium = TRUE requires column `crown_lad` to be numeric.", 
           call. = FALSE)
    }
    
    if (any(trees$crown_lad <= 0)) {
      stop(
        "turbid_medium = TRUE requires `crown_lad` to be strictly positive", 
        call. = FALSE)
    }
    
  } else {
    
    if (!"crown_openness" %in% names(trees)) {
      stop(
        "turbid_medium = FALSE requires column `crown_openness` in tree inventory.",
        call. = FALSE)
    }
    
    if (any(is.na(trees$crown_openness))) {
      stop(
        "turbid_medium = FALSE requires to define `crown_openness` for all trees.", 
        call. = FALSE)
    }
    
    if (!is.numeric(trees$crown_openness)) {
      stop("turbid_medium = FALSE requires column `crown_openness` to be numeric.", 
           call. = FALSE)
    }
    
    if (any(trees$crown_openness < 0 | trees$crown_openness > 1)) {
      stop(
        "turbid_medium = FALSE requires `crown_openness` to be in [0,1]", 
        call. = FALSE)
    }
  }
  
  invisible(TRUE)
}