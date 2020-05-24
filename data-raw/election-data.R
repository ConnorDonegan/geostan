library(readr)
library(stringr)
library(sf)
path <- "~/repo/rhs-esf-sim-study/article/revisions/data/turmoil.rds"
ohio <- read_rds(path)
ohio <- ohio[,-str_which(names(ohio), "^w\\.")]
names(ohio)[str_which(names(ohio), "density")] <- "pop_density"
ohio <- dplyr::select(ohio,
               GEOID, county, gop_growth, historic_gop, trump_2016, total_2016,
               population, pop_density, college_educated, white_nonhispanic, unemployment)


## custom albers equal area projection
albers <- "+proj=aea +lat_1=40.68 +lat_2=39.7 +lat_0=40.18 +lon_0=-82.64 +x_0=400000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ohio <- st_transform(ohio, crs = albers)

## save data with package
usethis::use_data(ohio, overwrite = TRUE)
