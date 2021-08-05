iter=10
silent = TRUE
source("helpers.R")

context("stan_car")

test_that("Poisson CAR model works", {
    data(sentencing)
  SW(     fit <- stan_car(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    car_parts = prep_car_data(shape2mat(sentencing, "B")),
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    silent = silent)
  )
    expect_geostan(fit)    
})

test_that("CAR accepts covariate ME, mixed (un-) bounded", {
    data(georgia)
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance + ICE,
                        data = georgia,
                        C = shape2mat(georgia),
                        ME = list(se = data.frame(insurance = georgia$insurance.se,
                                                  ICE = georgia$ICE.se),
                        bounded = c(1, 0),
                        bounds = c(0, 100)
                        ),        
                        chains = 1,
                        iter = iter,
                        silent = silent)
    )
    expect_geostan(fit)
})



test_that("CAR accepts covariate ME with WX, mixed ME-non-ME", {
    data(georgia)
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance + ICE,
                        slx = ~ insurance + ICE,
                        data = georgia,
                        C = shape2mat(georgia),
                        ME = list(se = data.frame(insurance = georgia$insurance.se),
                        bounded = 0,
                        bounds = c(0, 100)
                        ),        
                        chains = 1,
                        iter = iter,
                        silent = silent)
    )
    expect_geostan(fit)
})

