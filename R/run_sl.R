#' Run SamsaRaLight radiative balance
#'
#' Computes light interception and radiative balance for a forest stand using
#' the SamsaRaLight ray-tracing engine.
#'
#' This function is the standard user interface to SamsaRaLight. Advanced
#' ray-tracing, sky-discretization, and interception parameters are internally
#' set to standard values and do not need to be provided.
#'
#' @inheritParams run_sl_advanced
#'
#' @details
#' Internally, \code{run_sl()} calls \link{run_sl_advanced} using standard
#' ray-tracing, sky-discretization, and interception parameters:
#' \itemize{
#'   \item torus borders are enabled;
#'   \item crowns are represented as turbid media;
#'   \item the leaf extinction coefficient is \code{0.5};
#'   \item the clumping factor is \code{1};
#'   \item trunk interception is enabled;
#'   \item the minimum ray altitude angle is \code{10} degrees;
#'   \item direct rays use a starting offset of \code{0} degrees and an
#'     angular step of \code{5} degrees;
#'   \item diffuse rays use an angular step of \code{15} degrees;
#'   \item the Standard Overcast Sky model is used;
#'   \item the full year (days 1--365) is simulated.
#' }
#'
#' These defaults are intended for standard SamsaRaLight simulations. Users
#' requiring control over these parameters can use \link{run_sl_advanced}
#' directly.
#'
#' When \code{parallel_mode = TRUE}, SamsaRaLight uses OpenMP to parallelise
#' ray--target computations. BLAS is configured to use a single thread during
#' the simulation to avoid competing parallel execution.
#'
#' @return An object of class \code{"sl_output"} (a list) containing:
#' \itemize{
#'   \item \code{output}: A list containing the simulation results:
#'     \itemize{
#'       \item \code{light}: A list containing light-interception results for
#'         trees, ground cells, and sensors.
#'       \item \code{monthly_rays}: The monthly ray discretization and
#'         associated radiation energies. Returned only when
#'         \code{detailed_output = TRUE}.
#'       \item \code{interceptions}: Detailed interception matrices. Returned
#'         only when \code{detailed_output = TRUE}.
#'     }
#'   \item \code{params}: A list containing the simulation parameters,
#'     including the simulated period, sky model, ray-discretization
#'     parameters, and interception-model parameters.
#'   \item \code{input}: A list containing the original \code{sl_stand} and
#'     \code{monthly_radiations} objects.
#' }
#'
#' When \code{detailed_output = FALSE}, the main output tables contain:
#' \itemize{
#'   \item \code{light$sensors}: \code{id_sensor}, \code{e}, \code{pacl},
#'     and \code{punobs};
#'   \item \code{light$cells}: \code{id_cell}, \code{e}, \code{pacl},
#'     and \code{punobs};
#'   \item \code{light$trees}: \code{id_tree}, \code{epot}, \code{e},
#'     \code{lci}, \code{eunobs}, and \code{rci}.
#' }
#'
#' @seealso
#' \link{create_sl_stand}, \link{check_inventory}, \link{check_sensors},
#' \link{get_monthly_radiations}, \link{check_monthly_radiations},
#' \link{run_sl_advanced}
#'
#' @examples
#' data_prenovel <- SamsaRaLight::data_prenovel
#'
#' stand <- create_sl_stand(
#'   trees_inv = data_prenovel$trees,
#'   sensors = data_prenovel$sensors,
#'   cell_size = 10,
#'   latitude = data_prenovel$info$latitude,
#'   slope = data_prenovel$info$slope,
#'   aspect = data_prenovel$info$aspect,
#'   north2x = data_prenovel$info$north2x
#' )
#'
#' out <- run_sl(
#'   sl_stand = stand,
#'   monthly_radiations = data_prenovel$radiations
#' )
#'
#' str(out)
#'
#' @export
#' 
run_sl <- function(
    sl_stand,
    monthly_radiations,
    sensors_only = FALSE,
    detailed_output = FALSE,
    parallel_mode = FALSE,
    n_threads = NULL,
    verbose = TRUE
) {
  
  run_sl_advanced(
    sl_stand = sl_stand,
    monthly_radiations = monthly_radiations,
    sensors_only = sensors_only,
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
    detailed_output = detailed_output,
    parallel_mode = parallel_mode,
    n_threads = n_threads,
    verbose = verbose
  )
}