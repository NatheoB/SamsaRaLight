#' Create planar (x, y) coordinates from longitude / latitude
#' 
#' This function converts geographic coordinates (\code{lon}, \code{lat})
#' into planar coordinates (\code{x}, \code{y})
#' using an automatically selected UTM projection.
#'
#' @param df A data.frame containing geographic coordinates (\code{lon} and \code{lat}).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{df}}{The input data.frame with added \code{x} and \code{y} columns (meters, UTM).}
#'   \item{\code{epsg}}{EPSG code of the UTM projection used.}
#' }
#'
#' @details
#' The input data.frame must contain geographic coordinates (\code{lon} and \code{lat}).
#' Planar coordinates (\code{x}, \code{y}) are automatically computed using the
#' UTM zone inferred from the mean longitude and hemisphere inferred from the mean latitude.
#'
#' @importFrom sf st_as_sf st_transform st_coordinates
#' @export
create_xy_from_lonlat <- function(df) {
  
  ## ---- check coordinates ----------------------------------------------
  needs_conversion <- check_coordinates(df, verbose = F)
  
  if (!needs_conversion) {
    stop(
      "The data.frame already contains planar coordinates (`x`, `y`). ",
      "No conversion is needed."
    )
  }
  
  ## ---- compute EPSG ----------------------------------------------------
  epsg_utm <- get_utm_epsg(mean(df$lon, na.rm = TRUE),
                           mean(df$lat, na.rm = TRUE))
  
  ## ---- convert data ----------------------------------------------------
  df <- convert_lonlat(df, epsg_utm)
  
  ## ---- output ----------------------------------------------------------
  list(
    df   = df,
    epsg = epsg_utm
  )
}


#' Get EPSG code for UTM zone from longitude and latitude
#'
#' @param lon Mean longitude (decimal degrees)
#' @param lat Mean latitude (decimal degrees)
#'
#' @return EPSG code (numeric)
#'
#' @keywords internal
get_utm_epsg <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) {
    32600 + zone
  } else {
    32700 + zone
  }
}


#' Convert data.frame from lon/lat to UTM planar coordinates
#'
#' @param df A data.frame containing \code{lon} and \code{lat} columns
#' @param epsg_utm EPSG code of the UTM projection
#'
#' @return Data.frame with added \code{x} and \code{y} columns (meters)
#'
#' @importFrom sf st_as_sf st_transform st_coordinates
#' @keywords internal
convert_lonlat <- function(df, epsg_utm) {
  sf_pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  sf_pts_utm <- st_transform(sf_pts, crs = epsg_utm)
  coords <- st_coordinates(sf_pts_utm)
  df$x <- coords[, 1]
  df$y <- coords[, 2]
  df
}