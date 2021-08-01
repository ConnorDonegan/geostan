library(tidyverse)
library(sf)

## load county mortality data from Donegan, Chun, Griffith 2021
d.url <- "https://raw.githubusercontent.com/ConnorDonegan/survey-HBM/main/data/county-data.csv"
df <- read.csv(d.url)
ga <- df[str_which(df$county, ", Georgia"), ]
ga$GEOID <- as.character(ga$GEOID)

## load county geometry from US Census Bureau
sp.url <- "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_20m.zip"
tmp.dir <- tempfile(fileext=".zip")
download.file(sp.url, destfile = tmp.dir)
unzip(tmp.dir, exdir = "counties")
us <- st_read("counties", "cb_2018_us_county_20m")
sf.ga <- dplyr::filter(us, STATEFP == "13")

## join mortality data to geometry
sf.ga <- inner_join(sf.ga, ga, by = "GEOID")
sf.ga <- dplyr::select(sf.ga,
                       -c(STATEFP, COUNTYFP, COUNTYNS, AFFGEOID, county, LSAD)
                       )
georgia <- sf.ga

## make availabe to geostan users
usethis::use_data(georgia, overwrite = TRUE)
