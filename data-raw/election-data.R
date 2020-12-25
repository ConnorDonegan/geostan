#############################################################
## Update Ohio election data from Donegan et al. 2020 to
## include ACS standard errors on population characteristics.
## Standard errors on percentages are listed in ACS data
## profile tables (DP02, DP03, DP05) at data.census.gov
#############################################################

library(tidyverse)
library(sf)

## load data from Donegan et al 2020 paper
ohio.address <- "https://github.com/ConnorDonegan/ESF/raw/main/ohio.rda"
download.file(ohio.address, destfile = "ohio.rda")
load("ohio.rda")

## keep election data and BLS unemployment estimate
ohio <- dplyr::select(ohio,
               GEOID, county, gop_growth, historic_gop, trump_2016, total_2016, unemployment)

## custom albers equal area projection
albers <- "+proj=aea +lat_1=40.68 +lat_2=39.7 +lat_0=40.18 +lon_0=-82.64 +x_0=400000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ohio <- st_transform(ohio, crs = albers)

## DP02 educational attainment
dp02 <- read_csv("DP02-acs-2016-Ohio-extract.csv")
names(dp02) <- c("GEOID", "college_educated", "college.moe")
dp02 <- dp02 %>%
    mutate(
        college_educated.se = college.moe / 1.645,
        GEOID = str_extract(GEOID, ".{5}$")
    ) %>%
    select(-c(college.moe))

## DP03 unemployment rate (differs from BLS estimate)
dp03 <- read_csv("DP03-acs-2016-Ohio-extract.csv")
names(dp03) <- c("GEOID", "unemployment.acs", "unemployment.acs.moe")
dp03 <- dp03 %>%
    mutate(
        unemployment.acs.se = unemployment.acs.moe / 1.645,
        GEOID = str_extract(GEOID, ".{5}$")
    ) %>%
    select(-c(unemployment.acs.moe))

## DP05 race and ethnicity (total population, percent non-hispanic white)
dp05 <- read_csv("DP05-acs-2016-Ohio-extract.csv")
names(dp05) <- c("GEOID", "population", "white_nonhispanic", "white.moe")
dp05 <- dp05 %>%
    mutate(
        white_nonhispanic.se = white.moe / 1.645,
        GEOID = str_extract(GEOID, ".{5}$")
    ) %>%
    select(-c(white.moe))

## join
ohio <- ohio %>%
    inner_join(dp02, by = "GEOID") %>%
    inner_join(dp03, by = "GEOID") %>%
    inner_join(dp05, by = "GEOID")

## population density
area <- units::set_units(st_area(ohio), mi^2)
ohio$pop_density <- as.numeric(ohio$population / area)

## order nicely
ohio <- ohio %>%
    select(
        GEOID, county,
        gop_growth, historic_gop, trump_2016, total_2016,
        population, pop_density,
        white_nonhispanic, white_nonhispanic.se,
        college_educated, college_educated.se,
        unemployment, unemployment.acs, unemployment.acs.se
    ) 

## save data with package
usethis::use_data(ohio, overwrite = TRUE)
unlink("ohio.rda")
