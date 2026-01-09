#' Validate monthly radiation input
#'
#' @description
#' Checks that a monthly radiation data.frame is correctly formatted and
#' physically valid for use by the light interception and ray-tracing model.
#' The table must contain exactly 12 months of radiation data.
#'
#' @param x A data.frame with monthly radiation values, typically produced by
#'   \code{\link{get_monthly_radiations}}.
#'   
#' @param verbose Logical; if \code{TRUE}, informative messages are printed.
#'
#' @details
#' The input must contain the following columns:
#' \describe{
#'   \item{month}{Integer month number (1-12)}
#'   \item{Hrad}{Monthly global horizontal irradiation (MJ/m2)}
#'   \item{DGratio}{Diffuse-to-global radiation ratio (unitless, 0-1)}
#' }
#'
#' The function checks:
#' \describe{
#'   \item{1}{The object is a data.frame.}
#'   \item{2}{Required columns are present.}
#'   \item{3}{There are exactly 12 months.}
#'   \item{4}{Each month (1-12) is present exactly once.}
#'   \item{5}{Data are numeric and finite.}
#'   \item{6}{Hrad strictly positive.}
#'   \item{7}{DGratio between 0 and 1.}
#'   \item{8}{Months are in increasing order.}
#' }
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#'
#' @export
check_monthly_radiations <- function(x, verbose = TRUE) {
  
  ## ---- basic structure ------------------------------------------------------
  if (!is.data.frame(x)) {
    stop("Monthly radiation input must be a data.frame.", call. = FALSE)
  }
  
  if (nrow(x) == 0) {
    stop("Monthly radiation data.frame is empty.", call. = FALSE)
  }
  
  ## ---- required columns -----------------------------------------------------
  required <- c("month", "Hrad", "DGratio")
  missing <- setdiff(required, names(x))
  
  if (length(missing) > 0) {
    stop("Missing required column(s): ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }
  
  ## ---- types ---------------------------------------------------------------
  if (!is.numeric(x$month)) stop("`month` must be numeric.", call. = FALSE)
  if (!is.numeric(x$Hrad)) stop("`Hrad` must be numeric.", call. = FALSE)
  if (!is.numeric(x$DGratio)) stop("`DGratio` must be numeric.", call. = FALSE)
  
  ## ---- finite --------------------------------------------------------------
  if (any(!is.finite(x$month))) stop("`month` contains non-finite values.", call. = FALSE)
  if (any(!is.finite(x$Hrad))) stop("`Hrad` contains non-finite values.", call. = FALSE)
  if (any(!is.finite(x$DGratio))) stop("`DGratio` contains non-finite values.", call. = FALSE)
  
  ## ---- months --------------------------------------------------------------
  if (nrow(x) != 12) {
    stop("Monthly radiation table must contain exactly 12 rows (one per month).", call. = FALSE)
  }
  
  if (any(x$month %% 1 != 0)) {
    stop("`month` must contain integer values.", call. = FALSE)
  }
  
  if (!setequal(x$month, 1:12)) {
    stop("`month` must contain exactly the values between 1 and 12, each once.", call. = FALSE)
  }
  
  if (anyDuplicated(x$month)) {
    stop("`month` contains duplicated values.", call. = FALSE)
  }
  
  ## ---- ordering ------------------------------------------------------------
  if (!all(x$month == sort(x$month))) {
    stop("`month` must be in increasing order from 1 to 12.", call. = FALSE)
  }
  
  ## ---- physical validity ---------------------------------------------------
  if (any(x$Hrad < 0)) {
    stop("`Hrad` must be non-negative.", call. = FALSE)
  }
  
  if (any(x$DGratio < 0 | x$DGratio > 1)) {
    stop("`DGratio` must be between 0 and 1.", call. = FALSE)
  }
  
  ## ---- success --------------------------------------------------------------
  if (verbose) message("Radiation table successfully validated.")
  invisible(TRUE)
}
