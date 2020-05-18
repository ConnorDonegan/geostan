## get Ohio election data with ACS estimates and BLS unemployment estimates
pkgs <- c("rgdal", "geostan", "sf", "spdep", "tidyverse", "tidycensus", "readxl")
lapply(pkgs, require, character.only = TRUE)

# ACS variable search
# V <- load_variables(2016, "acs5", cache = TRUE)
# V %>% filter(name == 'B15003_022')

## Total population, college educated population over 25, non-hispanic whites
year = 2016
survey = "acs5"
options(tigris_use_cache = TRUE)
variables = c('B01003_001', 'B15003_001', 'B15003_022', 'B03002_003') # total population, total population over 25, college educated over 25, non-hispanic white

census <- tidycensus::get_acs(
  geography = "county",
  state = "OH",
  variables = variables,
  moe_level = 95,
  survey = survey,
  output = "wide",
  geometry = TRUE,
  keep_geo_vars = FALSE,
  cb = TRUE,
  #shift_geo = TRUE,
  year = year)

census <- census %>%
  transmute(GEOID = GEOID,
            population = B01003_001E,
            college_educated = 100 * (B15003_022E / B15003_001E),
            white_nonhispanic = 100 * (B03002_003E / population)) 

## County unemployment estimates, 2016
tmp <- tempfile(fileext = ".xlsx")
bls_path <- "https://www.bls.gov/lau/laucnty16.xlsx"
download.file(bls_path, tmp)
bls <- read_excel(tmp, skip = 5, col_names = FALSE)
names(bls) <- c("LAUS", "STATEFP", "COUNTYFP", "Name", "year", "skip", "labor_force", "employed", "unemployed", "unemployment_rate")

bls <- bls %>%
  dplyr::transmute(GEOID = paste0(STATEFP, COUNTYFP),
                   unemployment = unemployment_rate)

# MIT Election Data and Science Lab, 2018, "County Presidential Election Returns 2000-2016",
#   https://doi.org/10.7910/DVN/VOQCHQ, Harvard Dataverse, V1, UNF:6:ZaxsDvp6RsFitno8ZYlK7w== [fileUNF]
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ
tmpdf <- tempfile(fileext = ".csv")
download.file("https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/VOQCHQ/UM1EZF",
              destfile = tmpdf)

county_pres <- read_tsv(tmpdf) %>%
  dplyr::filter(state_po == "OH" &
                  party %in% c("democrat", "republican")) %>%
  mutate(GEOID = stringr::str_pad(FIPS, width = 5, side = 'left', pad = '0')) %>%
  group_by(GEOID, year) %>% 
  mutate(totalvotes = sum(candidatevotes, na.rm = TRUE)) %>%
  filter(party == "republican") %>%
  ungroup 

## historic two party gop vote share
past <- county_pres %>%
  filter(year != 2016) %>%
  group_by(GEOID) %>%
  summarise(historic_gop = 100 * sum(candidatevotes) / sum(totalvotes)) %>%
  ungroup

## summarize 2016 results and join to the past
tmp_ohio <- county_pres %>%
  filter(year == 2016) %>%
  mutate(trump_2016 = candidatevotes,
         total_2016 = totalvotes,
         trump_pct = 100 * candidatevotes / totalvotes
  ) %>% 
  inner_join(past, by = "GEOID") %>%
  transmute(state_po,
            county,
            GEOID,
            gop_growth = trump_pct - historic_gop,
            historic_gop,
            trump_2016,
            total_2016)

## join census to bls and ohio
ohio <- census %>%
  right_join(tmp_ohio, by = "GEOID") %>%
  left_join(bls, by = "GEOID") %>%
  dplyr::select(GEOID, state_po, county, gop_growth, historic_gop, trump_2016, total_2016, population, college_educated, white_nonhispanic, unemployment)

## population density
## see https://www.census.gov/quickfacts/fact/note/US/LND110210
get.shp <- function(url, folder = "shape") {
# provide url to download a shapefile and unzip it into the working directory
        tmp.dir <- tempfile(fileext=".zip")
	download.file(url, destfile = tmp.dir)
	unzip(tmp.dir, exdir = folder)
	list.files(folder, full.names = TRUE)
}
path = "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_5m.zip"
get.shp(path)
shp <- readOGR("shape", "cb_2018_us_county_5m")
shp <- shp@data[,c("GEOID", "ALAND")]
## mi^2
shp$ALAND <- as.numeric(shp$ALAND) / 2589988
ohio <- left_join(ohio, shp, by = "GEOID")
ohio$pop_density <- ohio$population / ohio$ALAND
ohio <- select(ohio,
               GEOID, ALAND, county, gop_growth, historic_gop, trump_2016, total_2016,
               population, pop_density, college_educated, white_nonhispanic, unemployment)


## custom albers equal area projection
albers <- "+proj=aea +lat_1=40.68 +lat_2=39.7 +lat_0=40.18 +lon_0=-82.64 +x_0=400000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ohio <- st_transform(ohio, crs = albers)

## save data with package
usethis::use_data(ohio, overwrite = TRUE)

## delete junk
unlink("shape", recursive = TRUE)
