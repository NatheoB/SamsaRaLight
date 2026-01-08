#' Compute direct and diffuse rays of a growing season
#'
#' @param monthly_rad data.frame - Monthly horizontal radiation (Hrad) and diffuse to global ratio (DGratio)
#'    Computed with function samsaRa::sl_get_monthlyrad()
#' @param latitude double - Latitude of the plot (in degrees)
#' @param start_day integer between 1 and 365 - First day of the vegetative period
#' @param end_day integer between 1 and 365 - Last day of the vegetative period
#' @param soc boolean - Standard Overcast Sky, if false: Uniform Overcast Sky
#' @param slope double - Slope of the plot (in degrees)
#' @param north_to_x_cw double - Angle from North to x axis clockwise. (in degrees)
#'    Default correspond to a Y axis oriented toward the North.
#' @param aspect double - Angle of slope bottom on the compass from the North, clockwise rotation (in degrees)
#'    northern aspect : 0, eastern aspect : 90, southern aspect : 180, western aspect : 270
#' @param height_anglemin double - Angle minimum between beam and soil (in degrees)
#' @param direct_startoffset double - Angle at which to start first direct ray (in degrees)
#' @param direct_anglestep double - Hour angle between two direct beams (in degrees)
#' @param diffuse_anglestep double - Hour angle between two diffuse beams (in degrees)
#'
#' @return list of 3 elements : horizontal energy (double), slope energy (double)
#' and rays (data.frame) with n rows and 5 columns:
#' \itemize{
#'  \item{azimut}{Azimut of the ray in radians}
#'  \item{height_angle}{Angle between beam and soil (in radians)}
#'  \item{e}{Energy of ray before crossing the canopy (in MJ.m-2)}
#'  \item{direct}{true if the ray is direct false if it is diffuse}
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom data.table fifelse
#'
#' @keywords internal
create_sl_rays <- function(monthly_rad,
                           latitude,
                           start_day = 1, end_day = 365,
                           soc = TRUE,
                           slope = 0,
                           north_to_x_cw = 90,
                           aspect = 0,
                           height_anglemin = 10,
                           direct_startoffset = 0,
                           direct_anglestep = 5,
                           diffuse_anglestep = 15) {
  
  ### Set global variables
  
  # Declination angle in degrees for each month
  DECLINATION_DEG <- c(-20.8, -12.7, -1.9, 9.9, 18.9, 23.1, 21.3, 13.7, 3.0, -8.8, -18.4, -23)
  
  
  ### Compute base variables ----
  latitude_rad            <- deg2rad(latitude)
  slope_rad               <- deg2rad(slope)
  north_to_x_cw_rad       <- deg2rad(north_to_x_cw)
  height_anglemin_rad     <- deg2rad(height_anglemin)
  direct_startoffset_rad  <- deg2rad(direct_startoffset)
  direct_anglestep_rad    <- deg2rad(direct_anglestep)
  diffuse_anglestep_rad   <- deg2rad(diffuse_anglestep)
  
  # Azimut of the vector orthogonal to the ground in the (x,y) system
  bottom_azimut_rad <- deg2rad(-aspect + north_to_x_cw)
  
  # Azimuth of south counterclockwise from x axis
  southazimut_ccw_rad <- pi + north_to_x_cw_rad
  southazimut_ccw_rad <- data.table::fifelse(southazimut_ccw_rad > 2 * pi,
                                             southazimut_ccw_rad - 2 * pi, southazimut_ccw_rad)
  
  # Declination angle in radians:
  # angle between the equator and a line drawn from the centre of the Earth to
  # the centre of the sun
  declination_rad <- deg2rad(DECLINATION_DEG)
  
  
  ### Get monthly radiations ----
  # Proportion of days for each month included within the vegetative period
  prop_months_veget <- prop_months(start_day, end_day)
  
  # Compute monthly global, direct and diffuse energies
  nrj_global_months <- monthly_rad$Hrad * prop_months_veget
  nrj_diffuse_months <- nrj_global_months * monthly_rad$DGratio
  nrj_direct_months <- nrj_global_months - nrj_diffuse_months
  
  
  ### Create diffuse and direct rays ----
  
  # Diffuse rays
  total_diffuse <- sum(nrj_diffuse_months)
  diffuse_rays <- sl_create_rays_diffuse(soc, total_diffuse,
                                         height_anglemin_rad, diffuse_anglestep_rad,
                                         slope_rad, bottom_azimut_rad)
  
  # Direct rays
  direct_rays <- sl_create_rays_direct(latitude_rad, declination_rad,
                                       nrj_direct_months,
                                       height_anglemin_rad, direct_anglestep_rad,
                                       slope_rad,
                                       bottom_azimut_rad, southazimut_ccw_rad,
                                       direct_startoffset_rad)
  
  # Create rays and return only rays with stricly positive energy (within growing season)
  rays <- dplyr::bind_rows(direct_rays$rays, diffuse_rays$rays)
  rays <- rays[rays$e_incident > 0,]
  rays <- dplyr::bind_cols(id_ray = 1:nrow(rays), rays)
  
  row.names(rays) <- NULL
  
  
  ### Return all rays and energies ----
  return(list("energies" = c("slope_direct" = direct_rays$e_slope,
                             "slope_diffuse" = diffuse_rays$e_slope,
                             "horizontal_direct" = direct_rays$e_horizontal,
                             "horizontal_diffuse" = diffuse_rays$e_horizontal),
              "rays" = rays))
}


#' Compute energy of diffuse ray
#'
#' @description Calculating SamsaraLight rays diffuse Energy in MJ/m2 of a plane
#' perpendicular to beam ray direction for a Standard Overcast Sky or
#' Uniform Overcast Sky are possible
#'
#' @param soc boolean - Standard Overcast Sky, if false: Uniform Overcast Sky
#' @param total_diffuse double - Total diffuse energy on a horizontal plan (in kWh.m-2)
#' @param height_angle_rad double - Angle between beam and soil (in radians)
#' @param diffuse_anglestep_rad double - Hour angle between two diffuse beams (in radians)
#'
#' @return Energy per square meter of a horizontal plan (in MJ.m-2)
#'
#' @export
#' @keywords internal
#' 
sl_compute_nrj_diffuse <- function(soc, total_diffuse,
                                   height_angle_rad,
                                   diffuse_anglestep_rad) {
  
  # Compute base variables
  meridian_nb <- 2 * pi / diffuse_anglestep_rad
  
  height_a_inf <- height_angle_rad - diffuse_anglestep_rad / 2
  sin_inf <- sin(height_a_inf)
  
  height_a_sup <- height_angle_rad + diffuse_anglestep_rad / 2
  sin_sup <- sin(height_a_sup)
  
  energy <- 0
  if (!soc) { # Uniform Overcast Sky, per square meter of a horizontal plan
    energy <- (2 * total_diffuse / meridian_nb) * (sin_sup^2 - sin_inf^2) / 2
    
  } else { # Standard Overcast Sky, Energy per square meter of a horizontal plan
    energy <- (6 * total_diffuse / (7 * meridian_nb)) *
      ((sin_sup^2 - sin_inf^2) / 2 + 2 * (sin_sup^3 - sin_inf^3) / 3)
  }
  energy <- energy / sin(height_angle_rad);
  
  return(energy)
}


#' Create diffuse rays
#'
#' @description Create SamsaraLight diffuse rays of a plane with a classical sky
#' hemisphere divided by meridians and parallels
#' can use Standard Overcast Sky or Uniform Overcast Sky
#'
#' @param soc boolean - Standard Overcast Sky, if false: Uniform Overcast Sky
#' @param total_diffuse double - Total diffuse energy on a horizontal plan (in kWh.m-2)
#' @param height_anglemin_rad double - Angle minimum between beam and soil (in radians)
#' @param diffuse_anglestep_rad double - Hour angle between two diffuse beams (in radians)
#' @param slope_rad double - Slope of the plan (in radians)
#' @param bottom_azimut_rad double - Azimut of the vector orthogonal to the ground
#' in the (x,y) system (in radians)
#'
#' @return list of 3 elements : horizontal energy (double), slope energy (double)
#' and diffuse rays (data.frame)
#'
#' @export
#' @keywords internal
#' 
sl_create_rays_diffuse <- function(soc, total_diffuse,
                                   height_anglemin_rad,
                                   diffuse_anglestep_rad,
                                   slope_rad,
                                   bottom_azimut_rad) {
  
  # All possible angles for diffuse rays (in radians)
  possible_angles_rad <- seq(from = diffuse_anglestep_rad / 2,
                             to = pi / 2,
                             by = diffuse_anglestep_rad)
  
  # All possible azimuts for diffuse rays (in radians)
  possible_azimuts_rad <- seq(from = diffuse_anglestep_rad / 2,
                              to = 2 * pi,
                              by = diffuse_anglestep_rad)
  
  # Compute energy of all possible rays (height_angles as a double vector)
  # BE CAREFUL : there can be round problems with transformation of angleStep
  # from degrees to radians and the last azimut can be very close to
  # 360 (one extra azimut)
  e_rays <- sl_compute_nrj_diffuse(soc, total_diffuse,
                                   possible_angles_rad, diffuse_anglestep_rad)
  
  # Create all combinations of possible angles of rays (in radians) and azimuts (in radians)
  ha_az <- data.frame(height_angle = possible_angles_rad,
                      e_incident = e_rays)
  ha_az <- ha_az[rep(seq_len(nrow(ha_az)), each = length(possible_azimuts_rad)), ]
  ha_az$azimut <- rep(possible_azimuts_rad, times = length(possible_angles_rad))
  
  # The Horizontal Diffuse reference is calculated for heightAngles > angleMin
  # For each azimut
  horizontal_angles <- ha_az$height_angle > height_anglemin_rad
  horizontal_diffuse <- sum( ha_az$e_incident[horizontal_angles] * sin( ha_az$height_angle[horizontal_angles] ) )
  
  # A beam is created only if it reaches the soil with an angle > angle_min
  # the cosinus of the angle between the vector orthogonal to slope and the beam
  # (given by scalar above) must be higher than sin(angle_min)
  ha_az$scalar <- cos(slope_rad) * sin(ha_az$height_angle) + sin(slope_rad) *
    cos(ha_az$height_angle) * cos(ha_az$azimut - bottom_azimut_rad);
  
  # Get beams to creates (remove beams with scalar <= sin(angle_min))
  ha_az <- ha_az[ha_az$scalar > sin(height_anglemin_rad),]
  
  # Compute slope diffuse energy
  slope_diffuse <- sum( ha_az$scalar * ha_az$e_incident )
  
  # Create the final diffuse rays dataframe
  ha_az$direct <- FALSE
  
  rays_diffuse <- ha_az[,c("azimut", "height_angle",
                           "e_incident", "direct")]
  
  # Return a list with rays and slope/horizontal energies
  return(list("e_horizontal" = horizontal_diffuse,
              "e_slope" = slope_diffuse,
              "rays" = rays_diffuse))
}


#' Create direct rays
#'
#' @description Create SamsaraLight diffuse rays of a plane with a classical sky
#' hemisphere divided by meridians and parallels
#' can use Standard Overcast Sky or Uniform Overcast Sky
#'
#' @param latitude_rad double - Latitude of the plot (in radians)
#' @param declination_rad double - Declination angle in radians:
#   angle between the equator and a line drawn from the centre of the Earth to
#   the centre of the sun
#' @param nrj_direct_months 12 double vect - Monthly direct energy on a horizontal plan (in kWh.m-2)
#' @param height_anglemin_rad double - Angle minimum between beam and soil (in radians)
#' @param direct_anglestep_rad double - Hour angle between two diffuse beams (in radians)
#' @param slope_rad double - Slope of the plan (in radians)
#' @param bottom_azimut_rad double - Azimut of the vector orthogonal to the ground
#' @param southazimut_ccw_rad double - Azimuth of south counterclockwise from x axis
#' in the (x,y) system (in radians)
#' @param direct_startoffset_rad double - Angle at which to start first direct ray (in radians)
#'
#' @return list of 3 elements : horizontal energy (double), slope energy (double)
#' and diffuse rays (data.frame)
#'
#' @export
#' @keywords internal
#' 
sl_create_rays_direct <- function(latitude_rad,
                                  declination_rad,
                                  nrj_direct_months,
                                  height_anglemin_rad,
                                  direct_anglestep_rad,
                                  slope_rad,
                                  bottom_azimut_rad,
                                  southazimut_ccw_rad,
                                  direct_startoffset_rad) {
  
  # Integrating sin(heightAngle) along the day
  sunrise_hourangle <- -acos(- tan(latitude_rad) * tan(declination_rad))
  sunset_hourangle <- - sunrise_hourangle
  
  day_sinheight_ang <- sin(latitude_rad) * sin(declination_rad) *
    (sunset_hourangle - sunrise_hourangle) + cos(latitude_rad) * cos(declination_rad) *
    (sin(sunset_hourangle) - sin(sunrise_hourangle))
  
  # Compute for each hour angle between [-pi, pi[
  hour_angles_rad <- seq(from = -pi + direct_startoffset_rad,
                         to = pi,
                         by = direct_anglestep_rad)
  hour_angles_rad <- hour_angles_rad[hour_angles_rad != pi]
  
  # For each combination of month (declination rad) and hour angles
  ha_dec <- data.frame(declination_rad = declination_rad,
                       energy_direct = nrj_direct_months,
                       day_sinheight_ang = day_sinheight_ang)
  
  ha_dec <- ha_dec[rep(seq_len(nrow(ha_dec)), each = length(hour_angles_rad)), ]
  ha_dec$hour_angle_rad <- rep(hour_angles_rad, times = length(declination_rad))
  
  # Compute height angle (in radians)
  ha_dec$height_angle_rad <- asin(sin(latitude_rad) * sin(ha_dec$declination_rad)
                                  + cos(latitude_rad) * cos(ha_dec$declination_rad) * cos(ha_dec$hour_angle_rad))
  
  # Compute sun azimut (in radians)
  ha_dec$azimut_rad <- sl_compute_sunazimut(latitude_rad,
                                            ha_dec$declination_rad,
                                            ha_dec$hour_angle_rad,
                                            ha_dec$height_angle_rad,
                                            southazimut_ccw_rad)
  
  # Compute hour sinus of height angle (in radians)
  ha_dec$hour_sin_height_ang <- sin(latitude_rad) * sin(ha_dec$declination_rad) * direct_anglestep_rad +
    cos(latitude_rad) * cos(ha_dec$declination_rad) *
    (sin(ha_dec$hour_angle + direct_anglestep_rad / 2) - sin(ha_dec$hour_angle - direct_anglestep_rad / 2))
  
  # Compute energy in MJ/m2 on a horizontal plane
  ha_dec$energy <- ha_dec$energy_direct * ha_dec$hour_sin_height_ang / ha_dec$day_sinheight_ang
  
  # The HorizontalDirect reference is calculated for heightAngles > angleMin
  horizontal_direct <- sum( ha_dec$energy[ha_dec$height_angle_rad > height_anglemin_rad] )
  
  # Compute perpendicular energy in MJ/m2 on a plane perpendicular to the ray
  # hour_sin_height_ang is very close to sin (height_angle)
  # so all rays have about the same amount of energy on the plane perpendicular to the ray.
  ha_dec$e_perpendicular <- ha_dec$energy / sin(ha_dec$height_angle_rad)
  
  # A ray is created only if it reaches the soil with an angle > angleMin
  # the cosinus of the angle between the vector orthogonal to slope and the ray
  # must be higher than sin(angleMin)
  # this cosinus is given by scalar :
  ha_dec$scalar <- cos(slope_rad) * sin(ha_dec$height_angle_rad) + sin(slope_rad) *
    cos(ha_dec$height_angle_rad) * cos(ha_dec$azimut_rad - bottom_azimut_rad)
  
  # Get beams to creates (remove beams with scalar <= sin(angle_min))
  ha_dec <- ha_dec[ha_dec$height_angle_rad > 0 & ha_dec$scalar > sin(height_anglemin_rad),]
  
  # Compute slope direct energy
  slope_direct <- sum( ha_dec$scalar * ha_dec$e_perpendicular )
  
  # Create the final diffuse rays dataframe
  ha_dec$e_incident <- ha_dec$e_perpendicular
  ha_dec$direct <- TRUE
  
  rays_direct <- ha_dec[,c("azimut_rad", "height_angle_rad",
                           "e_incident", "direct")]
  names(rays_direct) <- c("azimut", "height_angle", "e_incident", "direct")
  
  # Return a list with rays and slope/horizontal energies
  return(list("e_horizontal" = horizontal_direct,
              "e_slope" = slope_direct,
              "rays" = rays_direct))
}



#' Compute sun azimut
#'
#' @description Computaion of sun azimut for a given height angle reference system with
#' angle origin on X > 0 axis and trigonometric rotation
#'
#' @param latitude_rad double - Latitude of the plot (in radians)
#' @param declination_rad double - Declination angle in radians:
#'   angle between the equator and a line drawn from the centre of the Earth to
#'   the centre of the sun
#' @param hour_angle_rad double - Angle which the sun moves across the sky (in radians)
#' @param height_angle_rad double - Angle between beam and soil (in radians)
#' @param southazimut_ccw_rad double - Azimuth of south counterclockwise from x axis
#'  in the (x,y) system (in radians)
#'  
#' @importFrom data.table fifelse
#' 
#' @export
#' @keywords internal
#' 
sl_compute_sunazimut <- function(latitude_rad,
                                 declination_rad,
                                 hour_angle_rad,
                                 height_angle_rad,
                                 southazimut_ccw_rad) {
  
  # Solar position formulas in the reference system with azimut = 0 at
  # south and clockwise rotation
  sin_az <- cos(declination_rad) * sin(hour_angle_rad) / cos(height_angle_rad)
  
  cos_az <- (sin(latitude_rad) * cos(declination_rad) * cos(hour_angle_rad) -
               cos(latitude_rad) * sin(declination_rad)) / cos(height_angle_rad)
  
  azimut <- data.table::fifelse(cos_az >= 0,
                                asin(sin_az),
                                data.table::fifelse(sin_az >= 0,
                                                    pi - asin(sin_az),
                                                    -pi - asin(sin_az)))
  
  # Reference system with angle origin on X > 0 axis and trigonometric rotation
  # SouthAzimut_rad gives azimut of south direction in this system
  azimut <- southazimut_ccw_rad - azimut
  azimut <- azimut %% (2*pi)
  
  return(azimut)
}

#' Compute proportion of days for each month included between start and end days
#' @return 12-sized vector of double between 0 and 1 - Monthly proportion of growing days
#' @noRd
prop_months <- function(start_day, end_day) {
  
  # Number of days in each month (February always 28 days)
  NBMONTHDAYS <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  # Cumulative sum of the days across each month of the year
  cumsum_months <- cumsum(NBMONTHDAYS)
  
  # Number of days between start_day and last day of the month
  d1 <- cumsum_months - start_day + 1
  
  # Number of days between end_day and first day of the month
  d2 <- end_day - (cumsum_months - NBMONTHDAYS)
  
  # Compute proportion of days of each month included between start and end days
  abs(pmax(0, pmin(1,  1 - d1 / NBMONTHDAYS)) - pmax(0, pmin(1, d2 / NBMONTHDAYS)))
}