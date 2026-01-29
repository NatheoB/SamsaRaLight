#' Check coordinate columns and determine whether conversion is required
#'
#' This function checks whether a data frame contains valid planar
#' coordinates (`x`, `y`) or geographic coordinates (`lon`, `lat`).
#' If geographic coordinates are detected, the user is informed that
#' they must be converted to a planar coordinate system using
#' \code{create_xy_from_lonlat()}.
#'
#' @param df A data.frame containing spatial coordinates.
#' @param verbose Logical. If TRUE, informative messages are printed.
#'
#' @return Logical. Invisibly returns \code{TRUE} if coordinates must be converted
#'   from longitude/latitude to planar coordinates, and \code{FALSE}
#'   if planar coordinates are already present.
#'
#' @examples
#' df_xy <- data.frame(x = 1:3, y = 4:6)
#' check_coordinates(df_xy)
#'
#' df_lonlat <- data.frame(lon = c(10, 11), lat = c(45, 46))
#' check_coordinates(df_lonlat, verbose = TRUE)
#'
#' @export
check_coordinates <- function(df, verbose = TRUE) {
  
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }
  
  has_xy  <- all(c("x", "y") %in% names(df))
  has_lonlat  <- all(c("lon", "lat") %in% names(df))
  
  if (has_xy) {
    if (!is.numeric(df$x) || !is.numeric(df$y)) {
      stop("`x` and `y` columns must be numeric.", call. = FALSE)
    }
    
    if (verbose) message("Planar coordinates successfully validated.")
    
    return(invisible(FALSE))
  }
  
  if (has_lonlat) {
    if (!is.numeric(df$lon) || !is.numeric(df$lat)) {
      stop("`lon` and `lat` columns must be numeric.", call. = FALSE)
    }
    
    if (verbose) message("Geographic coordinates (`lon`, `lat`) detected.")
    
    return(invisible(TRUE))
  }
  
  stop(
    "Data must contain either planar coordinates (`x`, `y`) ",
    "or geographic coordinates (`lon`, `lat`).",
    call. = FALSE
  )
}
