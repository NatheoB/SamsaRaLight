# Create the SamsaraLight

Fetch monthly radiation data from PVGIS website (by API) between start
and end year (limit years are from 2005 to 2020). Fetched variables are
Hrad = horizontal plane irradiation and DGratio = ratio of diffuse to
global radiation (in horizontal plane).

! YOU NEED AN INTERNET CONNECTION TO ACCESS THE DATA BY API !

## Usage

``` r
get_monthly_radiations(latitude, longitude, start_year = 2005, end_year = 2020)
```

## Source

https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system_en

## Arguments

- latitude:

  latitude of the plot

- longitude:

  longitude of the plot

- start_year:

  positive integer between 2005 and 2020 - start year on which to fetch
  monthly data

- end_year:

  positive integer between 2005 and 2020 - end year on which to fetch
  monthly data

## Value

Monthly horizontal radiation (Hrad) and diffuse to global ratio
(DGratio) averaged between start_year and end_year
