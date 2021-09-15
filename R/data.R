#' Georgia all-cause, sex-specific mortality, ages 55-64, years 2014-2018
#'
#' @source
#'
#' Centers for Disease Control and Prevention, National Center for Health Statistics. Underlying Cause of Death 1999-2018 on CDC Wonder Online Database. 2020. Available online: \url{http://wonder.cdc.gov} (accessed on 19 October 2020).
#' 
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). ``Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure.'' *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Kyle Walker and Matt Herman (2020). tidycensus: Load US Census Boundary and Attribute Data as 'tidyverse' and 'sf'-Ready Data  Frames. R package version 0.11. \url{https://CRAN.R-project.org/package=tidycensus}
#'
#' US Census Bureau. Variance Replicate Tables, 2018. Available online: \url{https://www.census.gov/programs-surveys/acs/data/variance-tables.2018.html} (accessed on 19 October 2020).
#' 
#' @description A simple features (sf) object for Georgia counties with sex- and age-specific deaths and populations at risk (2014-2018), plus select estimates (with standard errors) of county characteristics. Standard errors of the ICE were calculated using the Census Bureau's variance replicate tables.
#' 
#' @format A simple features object with county geometries and the following columns:
#' \describe{
#'  \item{GEOID}{Six digit combined state and county FIPS code}
#'  \item{NAME}{County name}
#'  \item{ALAND}{Land area}
#'  \item{AWATER}{Water area}
#'  \item{population}{Census Bureau 2018 county population estimate}
#'  \item{white}{Percent White, ACS 2018 five-year estimate}
#'  \item{black}{Percent Black, ACS 2018 five-year estimate}
#'  \item{hisp}{Percent Hispanic/Latino, ACS 2018 five-year estimate}
#'  \item{ai}{Percent American Indian, ACS 2018 five-year estimate}
#'  \item{deaths.male}{Male deaths, 55-64 yo, 2014-2018}
#'  \item{pop.at.risk.male}{Population estimate, males, 55-64 yo, years 2014-2018 (total), ACS 2018 five-year estimate}
#'  \item{pop.at.risk.male.se}{Standard error of the pop.at.risk.male estimate}
#'  \item{deaths.female}{Female deaths, 55-64 yo, 2014-2018}
#'  \item{pop.at.risk.female}{Population estimate, females, 55-64 yo, years 2014-2018 (total), ACS 2018 five-year estimate}
#'  \item{pop.at.risk.female.se}{Standard error of the pop.at.risk.female estimate}
#'  \item{ICE}{Index of Concentration at the Extremes}
#'  \item{ICE.se}{Standard error of the ICE estimate, calculated using variance replicate tables}
#'  \item{income}{Median household income, ACS 2018 five-year estimate}
#'  \item{income.se}{Standard error of the income estimate}
#'  \item{college}{Percent of the population age 25 or higher than has a bachelors degree of higher, ACS 2018 five-year estimate}
#'  \item{college.se}{Standard error of the college estimate}
#'  \item{insurance}{Percent of the population with health insurance coverage, ACS 2018 five-year estimate}
#'  \item{insurance.se}{Standard error of the insurance estimate}
#'  \item{rate.male}{Raw (crude) age-specific male mortality rate, 2014-2018}
#'  \item{rate.female}{Raw (crude) age-specific female mortality rate, 2014-2018}
#'  \item{geometry}{simple features geometry for county boundaries}  
#' }
#' @examples
#' \dontrun{
#' library(sf)
#' data(georgia)
#' plot(georgia[,'rate.female'])
#' }
"georgia"


#' Florida state prison sentencing counts by county, 1905-1910
#'
#' @source Donegan, Connor. "The Making of Florida's 'Criminal Class': Race, Modernity and the Convict Leasing Program." Florida Historical Quarterly 97.4 (2019): 408-434. \url{https://osf.io/2wj7s/}.
#'
#' Mullen, Lincoln A. and Bratt, Jordon. "USABoundaries: Historical and Contemporary Boundaries of the United States of America,"
#'  Journal of Open Source Software 3, no. 23 (2018): 314, \doi{10.21105/joss.00314}.
#'
#' @description A spatial polygons data frame of historical 1910 county boundaries of Florida with aggregated state prison sentencing counts and census data.
#'  Sentencing and population counts are aggregates over the period 1905-1910, where populations were interpolated linearly between decennial censuses of 1900 and 1910.
#' @format A spatial polygons data frame with the following attributes:
#' \describe{
#'  \item{name}{County name}
#'  \item{wpop}{White population total for years 1905-1910}
#'  \item{bpop}{Black population total for years 1905-1910}
#'  \item{sents}{Number of state prison sentences, 1905-1910}
#'  \item{plantation_belt}{Binary indicator for inclusion in the plantation belt}
#'  \item{pct_ag_1910}{Percent of land area in agriculture, 1910}
#'  \item{expected_sents}{Expected sentences given demographic information and state level sentencing rates by race}
#'  \item{sir_raw}{Standardized incident ratio (observed/expected sentences)}
#' }
#' @examples
#' \dontrun{
#' data(sentencing)
#' }
"sentencing"

