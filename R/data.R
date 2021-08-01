#' Ohio Presidential election results and county characteristics
#'
#' @source
#'
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. Spatial Statistics. \url{https://doi.org/10.1016/j.spasta.2020.100450}
#' 
#' MIT Election Data and Science Lab, 2018, "County Presidential Election Returns 2000-2016", \url{https://doi.org/10.7910/DVN/VOQCHQ}, Harvard Dataverse, V1, UNF:6:ZaxsDvp6RsFitno8ZYlK7w== [fileUNF]
#' 
#' US Bureau of Labor Statistics, Local Force By County, 2016 Annual Averages. \url{https://www.bls.gov/lau/#cntyaa}
#'
#' US Census Bureau, 2016. American Community Survey. Tables DP02, DP03, DP05 5-year estimates.
#'
#' Kyle Walker, 2018. Tidycensus. \url{https://walkerke.github.io/tidycensus/index.html}
#'
#' @description A simple features spatial dataset containing county attributes and Presidential election data for Ohio.
#'  This is a processed version of the MIT "County Presidential Election Returns 2000-2016" data set augmented with data from the American Community Survey and Bureau of Labor Statistics.
#'
#' @details The Bureau of Labor Statistics unemployment estimates are based on the ACS estimates as well as additional information. The BLS estimates were used for the analysis in Donegan et al. (2020). The ACS estimates have standard errors (the BLS estimates do not).
#'
#' @format A simple feature collection including the following attributes:
#' \describe{
#'  \item{GEOID}{Six digit combined state and county FIPS code}
#'  \item{county}{County name}
#'  \item{gop_growth}{Change in the Republican vote share from historic (2000 - 2012) average vote share to 2016 vote share (i.e. (trump_2016/total_2016) - historic_gop). Data only includes Democratic and Republican votes (i.e. two-party vote share).}
#'  \item{historic_gop}{Average Republican share of all major party votes, 2000 to 2012}
#'  \item{trump_2016}{Number of votes for Donald Trump in 2016}
#'  \item{total_2016}{Total number of Democratic and Republican votes in 2016}
#'  \item{population}{ACS 2016 5 year population estimate}
#' \item{pop_density}{Population per square mile (population / ALAND)}
#'  \item{college_educated}{ACS 2016 5 year estimate of the percent of the population 25 and older with a bachelors degree or higher}
#' \item{college_educated.se}{Standard error for \code{college_educated}}
#'  \item{white_nonhispanic}{ACS 2016 5 year estimate of the percent of non-hispanic whites in the population}
#' \item{white_nonhispanic.se}{Standard error for \code{white_nonhispanic}}
#'  \item{unemployment}{Local area unemployment estimate from the Bureau of Labor Statistics, 2016}
#' \item{unemployment.acs}{Civilian unemployment rate estimate from the American Community Survey, 2016}
#' \item{unemployment.acs.se}{Standard error for \code{unemployment.acs}}
#' \item{geometry}{County boundaries (polygons) in simple features format}
#' }
#' @examples
#' \dontrun{
#' library(sf)
#' data(ohio)
#' }
"ohio"

#' Florida state prison sentencing counts by county, 1905-1910
#'
#' @source Donegan, Connor. "The Making of Florida's 'Criminal Class': Race, Modernity and the Convict Leasing Program." Florida Historical Quarterly 97.4 (2019): 408-434.
#'  \url{https://osf.io/2wj7s/}.
#'
#' Mullen, Lincoln A. and Bratt, Jordon. "USABoundaries: Historical and Contemporary Boundaries of the United States of America,"
#'  Journal of Open Source Software 3, no. 23 (2018): 314, \url{https://doi.org/10.21105/joss.00314}.
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


#' Georgia all-cause, sex-specific mortality, ages 55-64, years 2014-2018
#'
#' @source
#'
#' Centers for Disease Control and Prevention, National Center for Health Statistics. Underlying Cause of Death 1999-2018 on CDC Wonder Online Database. 2020. Available online: \url{http://wonder.cdc.gov/ucd-iid10.html} (accessed on 19 October 2020).
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
#'  \item{pop.at.risk.male}{Population estimate, males, 55-64 yo, years 2014-2018 (total), ACS 2018 five-year estimtae}
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
