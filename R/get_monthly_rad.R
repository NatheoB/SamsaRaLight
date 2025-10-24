#' Fetch monthly radiations
#'
#' @description Fetch monthly radiation data from PVGIS website (by API)
#'    between start and end year (limit years are from 2005 to 2020).
#'    Fetched variables are Hrad = horizontal plane irradiation and
#'    DGratio = ratio of diffuse to global radiation (in horizontal plane).
#'
#'    ! YOU NEED AN INTERNET CONNECTION TO ACCESS THE DATA BY API !
#'
#'
#' @param latitude latitude of the plot
#' @param longitude longitude of the plot
#' @param start_year positive integer between 2005 and 2020 - start year on which to fetch monthly data
#' @param end_year positive integer between 2005 and 2020 - end year on which to fetch monthly data
#' @param average_years if TRUE, years are averaged to obtain a single value for each month
#'
#' @return Monthly horizontal radiation (Hrad) and diffuse to global ratio (DGratio)
#'    averaged between start_year and end_year
#'
#' @source https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system_en
#'
#' @import readr httr
#' @importFrom dplyr recode mutate select arrange group_by summarize_all %>% 
#'
#' @export
get_monthly_rad <- function(latitude, longitude,
                            start_year = 2005,
                            end_year = 2020,
                            average_years = TRUE) {

  # Check for start and end years (PVGIS go from 2005 to 2020) ----
  if (start_year < 2005) stop("PVGIS data start from year 2005")
  if (end_year > 2020) stop("PVGIS data end at year 2020")
  if (start_year > end_year) stop("start_year must be lower or equal to end_year")


  # Get the data for each month and year, as string format ----
  out_str <- tryCatch(expr = {
    
    # Make our request to the API
    res <- httr::GET(paste0("https://re.jrc.ec.europa.eu/api/v5_2/MRcalc?", # MRcalc = Monthly radiation tool
                            "lat=", latitude, # lat latitude
                            "&lon=", longitude, # lon = longitude
                            "&startyear=", start_year, # Start year on which to fetch monthly data
                            "&endyear=", end_year, # End year on which to fetch monthly data
                            "&horirrad=1", # horirrad : Output horizontal plane irradiation
                            "&d2g=1", # d2g : Output monthly values of the ratio of diffuse to global radiation (horizontal plane)
                            "&outputformat=basic")) # outputformat as "basic" = output without text as csv
    
    # Transform request in binary into string
    rawToChar(res$content)
    
  }, error = function(e) {
    print(e)
    return(NULL)
  })
  
  # Check if problems
  if (grepl("message", out_str)) {
    print(out_str)
    stop("Wrong coordinates")
  }

  # Convert list string output into a single data.frame (binded by row)
  out_df <- read.table(text=out_str,
                       col.names = c("year", "month", "Hrad", "DGratio")) %>%
    
    # Rename months with integers between 1 and 12, and arrange from January to December
    dplyr::mutate(month = as.integer(
      dplyr::recode(month,
                    "Jan" = 1, "Feb" = 2, "Mar" = 3, "Apr" = 4,
                    "May" = 5, "Jun" = 6, "Jul" = 7, "Aug" = 8,
                    "Sep" = 9, "Oct" = 10, "Nov" = 11, "Dec" = 12))) %>%
    dplyr::mutate(Hrad = Hrad * 1000 * 3600 / 1e6) %>%  # Convert kWh into MJ
    dplyr::select(year, month, Hrad, DGratio) %>%
    dplyr::arrange(year, month)

  
  # Compute monthly mean of Hrad and DGratio between start and end years if specified
  if (average_years) {
    out_df <- out_df %>%
      dplyr::select(-year) %>%
      dplyr::group_by(month) %>%
      dplyr::summarize_all(mean)
  }

  # Convert to data.frame
  out_df <- as.data.frame(out_df)
  
  return(out_df)
}

