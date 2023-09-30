# download sf file with Florida counties and state prison sentencing data, 1905-1910
library(readr)
library(sf)
tf <- tempfile(fileext = ".rds")
download.file("https://github.com/ConnorDonegan/convict-leasing/raw/master/data/constructed/florida-1910-sp.rds", tf)
fl <- read_rds(tf)
fl@data <- fl@data[, c("name", "wpop", "bpop", "sents", "plantation_belt", "pct_ag_1910", "expected_sents")]
fl@data$sir_raw <- with(fl@data, sents / expected_sents)
sentencing <- st_as_sf(fl)
usethis::use_data(sentencing, overwrite = TRUE)

