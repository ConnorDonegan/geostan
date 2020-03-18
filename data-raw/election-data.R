library(rgdal)
library(sf)
library(spdep)
library(tidyverse)
library(tidycensus)
library(readxl)

# ACS variable search
# V <- load_variables(2016, "acs5", cache = TRUE)
# V %>% filter(name == 'B15003_022')

# Total population, college educated population over 25, non-hispanic whites
year = 2016
survey = "acs5"
options(tigris_use_cache = TRUE)
variables = c('B01003_001', 'B15003_001', 'B15003_022', 'B03002_003') # total population, total population over 25, college educated over 25, non-hispanic white

census <- tidycensus::get_acs(
  geography = "county",
  state = c("WI", "IL", "IN", "OH", "MI", "PA"),
  variables = variables,
  moe_level = 95,
  survey = survey,
  output = "wide",
  geometry = TRUE,
  keep_geo_vars = FALSE,
  cb = TRUE,
  #shift_geo = TRUE,
  year = year)

# don't like the hanlding of water in this shape file. Will us county boundary file instead for geometry
census <- census %>%
  transmute(GEOID = GEOID,
            population = B01003_001E,
            college_educated = 100 * (B15003_022E / B15003_001E),
            white_nonhispanic = 100 * (B03002_003E / population)) %>%
  mutate(population = population / 1e3)

# County unemployment estimates, 2016
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
  dplyr::filter(state_po %in% c("WI", "IL", "IN", "OH", "MI", "PA")) %>%
  mutate(GEOID = stringr::str_pad(FIPS, width = 5, side = 'left', pad = '0')) %>%
  group_by(GEOID, year) %>%
  mutate(totalvotes = sum(candidatevotes, na.rm = TRUE)) %>%
  ungroup

past <- county_pres %>%
  filter(party == "republican" &
           year != 2016 &
           !is.na(state_po)) %>%
  group_by(GEOID) %>%
  summarise(historic_gop = 100 * sum(candidatevotes, na.rm = TRUE) / sum(totalvotes, na.rm = TRUE)) %>%
  ungroup()

tmp_turmoil <- county_pres %>%
  filter(year == 2016 &
           party == "republican" &
           !is.na(state_po)) %>%
  mutate(trump_2016 = candidatevotes,
         total_2016 = totalvotes,
         trump_pct = 100 * candidatevotes / totalvotes
         ) %>%
  inner_join(past, by = "GEOID") %>%
  transmute(state_po,
            county,
            GEOID, #= recode(GEOID, `46113` = "46102"),
            gop_growth = trump_pct - historic_gop,
            historic_gop,
            trump_2016,
            total_2016)

# join together
turmoil <- census %>%
  right_join(tmp_turmoil, by = "GEOID") %>%
  left_join(bls, by = "GEOID") %>%
  dplyr::select(GEOID, state_po, county, gop_growth, historic_gop, trump_2016, total_2016, population, college_educated, white_nonhispanic, unemployment)

turmoil <- st_transform(turmoil, 3174)
usethis::use_data(turmoil, overwrite = TRUE)

