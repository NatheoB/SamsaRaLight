#' Run SamsaRaLight radiative balance
#'
#' This function computes light interception and radiative balance for a forest
#' stand using the SamsaRaLight ray-tracing engine.
#'
#' It is the **standard user interface** of SamsaRaLight.  
#' Advanced ray-tracing and sky discretization parameters are internally set to
#' robust defaults and do not need to be provided.
#'
#' @param sl_stand An object of class \code{"sl_stand"} describing the forest stand,
#'   created with \link{create_sl_stand}. It contains trees, sensors, terrain,
#'   and grid geometry.
#'
#' @param monthly_radiations A data.frame of monthly horizontal radiation
#'   (\code{Hrad}, in MJ m\eqn{^{-2}}) and diffuse-to-global ratio
#'   (\code{DGratio}), typically obtained using
#'   \link{get_monthly_radiations} and checked using \link{check_monthly_radiations}..
#'
#' @param latitude Numeric. Latitude of the stand (degrees).
#'
#' @param sensors_only Logical.  
#' If \code{TRUE}, compute light interception only for sensors (much faster).
#'
#' @param use_torus Logical.  
#' If \code{TRUE}, stand borders are treated using a torus (periodic) geometry,
#' mimicking an infinite forest, representative of the virtual stand.  
#' If \code{FALSE}, borders are open (surrounded by grassland).
#'
#' @param turbid_medium Logical.  
#' If \code{TRUE}, tree crowns are treated as a turbid medium (Beer–Lambert law
#' using \code{crown_lad}).  
#' If \code{FALSE}, crowns are treated as porous envelopes (using \code{crown_openness}).
#'
#' @param detailed_output Logical.  
#' If \code{TRUE}, the output contains detailed diffuse/direct energies in the \code{light} datasets, full
#' interception matrices \code{interceptions} and output of ray discretization \code{monthy_rays}.  
#' If \code{FALSE}, only total energies are returned (recommended for most uses).
#'
#' @param parallel_mode logical. If TRUE, ray–target computations are parallelised
#'   using OpenMP. If FALSE, the model runs in single-thread mode. SamsaRaLight uses OpenMP for ray–target parallelisation. To avoid competition
#'   between OpenMP and BLAS (matrix algebra libraries), BLAS is automatically forced
#'   to single-thread mode during the simulation. Using \code{parallel_mode = TRUE} is strongly recommended for large stands
#'   or fine ray discretisation, as computation time scales almost linearly with
#'   the number of available CPU cores.
#'   
#' @param n_threads integer or NULL. Number of CPU threads to use when
#'   \code{parallel_mode = TRUE}. If NULL (default), OpenMP automatically selects
#'   the number of available cores. If provided, must be a positive integer.
#'
#' @details
#' Internally, \code{run_sl()} calls the advanced engine
#' \code{run_sl_advanced()} with fixed ray-tracing and sky discretization.
#'
#' You should normally **not** use \code{SamsaRaLight:::run_sl_advanced()} directly unless you
#' are developing new ray-tracing configurations or doing methodological work.
#'
#' @return An object of class \code{"sl_output"}, containing:
#' \describe{
#'   \item{light}{
#'     A list of data.frames with simulated light interception:
#'     \itemize{
#'       \item \code{trees}: light intercepted by trees
#'       \item \code{cells}: light received by ground cells
#'       \item \code{sensors}: light received by sensors
#'     }
#'   }
#'   \item{info}{
#'     A list of metadata about the simulation (latitude, sky type, torus use, etc.).
#'   }
#'   \item{monthly_rays}{(only if \code{detailed_output = TRUE}) Discretization of monthly radiations}
#'   \item{interceptions}{(only if \code{detailed_output = TRUE}) interception matrices between trees and rays for each cell/sensor}
#' }
#'
#' @seealso
#' \link{create_sl_stand}, \link{check_inventory}, \link{check_sensors},
#' \link{get_monthly_radiations}, \link{check_monthly_radiations}
#'
#' @examples
#' \dontrun{
#' data_prenovel <- SamsaRaLight::data_prenovel
#'
#' stand <- create_sl_stand(
#'   trees = data_prenovel$trees,
#'   sensors = data_prenovel$sensors,
#'   cell_size = 5,
#'   slope = 10,
#'   aspect = 180,
#'   north2x = 0
#' )
#'
#' rad <- get_monthly_radiations(45.8, 3.1)
#'
#' out <- run_sl(
#'   sl_stand = stand,
#'   monthly_radiations = rad,
#'   latitude = 45.8
#' )
#'
#' out$light$trees
#' }
#'
#' @export
run_sl <- function(
    sl_stand,
    monthly_radiations,
    latitude,
    sensors_only = FALSE,
    use_torus = TRUE,
    turbid_medium = TRUE,
    detailed_output = FALSE,
    parallel_mode = FALSE,
    n_threads = NULL
) {
  
  SamsaRaLight:::run_sl_advanced(
    sl_stand = sl_stand,
    monthly_radiations = monthly_radiations,
    latitude = latitude,
    sensors_only = sensors_only,
    use_torus = use_torus,
    turbid_medium = turbid_medium,
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
    n_threads = n_threads
  )
}