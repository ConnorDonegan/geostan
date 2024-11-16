##
##
## A Monte Carlo study to test geostan implementation of spatial lag and error models
##
## Estimates should have following mathematical properties:
##   Non-biased
##   Monte Carlo rmse = reported SE 
##   Monte Carlo interval coverage rates = nominal coverage rates
##

#devtools::load_all("~/dev/geostan")
library(geostan)

## no. iterations
M = 15

# generate y using SLM or SEM
#sar_type <- commandArgs(trailingOnly = TRUE)
#sar_type <- match.arg(sar_type, c("SEM", "SLM"))
sar_type <- "SEM"

## a regular grid
ncol = 20
nrow = 20
sars <- prep_sar_data2(ncol, nrow, quiet = TRUE)
cars <- prep_car_data2(ncol, nrow, quiet = TRUE)
W <- sars$W
N <- nrow(W)

## View the grid by uncommenting and running:
#library(sf)
#sfc = sf::st_sfc(st_polygon(list(rbind(c(0,0), c(ncol,0), c(ncol,nrow), c(0,0)))))
#geom <- sf::st_make_grid(sfc, cellsize = 1, square = TRUE)
#df <- data.frame(geometry = geom)
#grid <- sf::st_as_sf(df)
# plot(grid)

# iterate:
b = 0.5
g = -0.5
rho.x <- 0.6
sigma.y = 0.3
rho.y = 0.6
pars <- c(const = 0, gamma = g, beta = b, rho = rho.y, scale = sigma.y)

res <- sapply(1:M, FUN = function(i) {
    
    x <- geostan::sim_sar(rho = rho.x, w = W, type = "SEM", quick = TRUE)
    wx <- as.numeric( W  %*% x )
    y <- geostan::sim_sar(mu = g * wx + b * x, rho = rho.y, w = W, sigma = sigma.y, type = sar_type, quick = TRUE)
    dat <- cbind(x = x , y = y)   

    fit_sem <- geostan::stan_sar(y ~ x,
                                 data = dat,
                                 type = "SDEM",
                                 sar_parts = sars,
                                 chains = 1,
                                 iter = 500,
                                 quiet = TRUE,
                                 slim = TRUE
                                 ) |>
        suppressWarnings()

    fit_car <- geostan::stan_car(y ~ x,
                                 slx = ~ x,
                                 data = dat,
                                 car_parts = cars,
                                 chains = 1,
                                 iter = 500,
                                 quiet = TRUE,
                                 slim = TRUE
                                 ) |>
        suppressWarnings()

    fit_glm <- geostan::stan_glm(y ~ x,
                                 slx = ~ x,
                                 data = dat,
                                 C = cars$C,
                                 chains = 1,
                                 iter = 500,
                                 quiet = TRUE,
                                 slim = TRUE
                                 ) |>
        suppressWarnings()
    
    x <- c(SEM = fit_sem$summary[c('intercept', 'w.x', 'x', 'sar_rho', 'sar_scale'), 'mean'],
    CAR = fit_car$summary[c('intercept', 'w.x', 'x', 'car_rho', 'car_scale'), 'mean'],
    GLM = fit_glm$summary[c('intercept', 'w.x', 'x', 'rho', 'sigma'), 'mean'])
    
    return (x)
    
    })

cat("\n**\nEstimates (rounded) should be close to SEM-DGP parameters, with allowance for CAR/GLM rho and scale\n",M, "iterations\n**\n")

est <- apply(res, 1, mean) |>
               round(1)

out <- data.frame(
    Mod = attributes(est)$names,
    Par = attributes(pars)$names,
    DGP = pars,
    Est = est,
    Error = est - pars
)

out |>
    print() |>
    suppressWarnings()

write.csv(out, "check-SEM-CAR-GLM-estimates-monte-carlo-output.csv")
