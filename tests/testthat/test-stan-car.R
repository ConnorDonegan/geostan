iter=10
refresh = 0
source("helpers.R")

context("stan_car")

test_that("Poisson CAR model works", {
    data(sentencing)
    SW(
        fit <- stan_car(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    car_parts = prep_car_data(shape2mat(sentencing, "B")),
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
  )
    expect_geostan(fit)    
})

test_that("CAR accepts covariate ME", {
    data(georgia)
    SW(
        fit <- stan_car(log(rate.male) ~ insurance + ICE,
                        data = georgia,
                        ME = prep_me_data(se = data.frame(insurance = georgia$insurance.se,
                                                          ICE = georgia$ICE.se)
                        ),
                        car_parts = prep_car_data(shape2mat(georgia, "B")),                        
                        chains = 1,
                        iter = iter,
                        refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("CAR accepts covariate ME with logit transform", {
    data(georgia)
    georgia$income <- georgia$income/1e3
    georgia$income.se <- georgia$income.se/1e3
    georgia$log_income <- log(georgia$income)
    georgia$log_income.se <- se_log(georgia$income, georgia$income.se)
    georgia$college <- georgia$college/1e3
    georgia$college.se <- georgia$college.se/1e3
    
    ME <- prep_me_data(se = data.frame(college = georgia$college.se,
                                       log_income = georgia$log_income.se),
                       logit = c(TRUE, FALSE),
                       bounds =c (0, Inf)
                       )
    SW(
        fit <- stan_car(log(rate.male) ~ college + log_income,
                        data = georgia,
                        ME = ME,
                        car_parts = prep_car_data(shape2mat(georgia, "B")),                        
                        chains = 1,
                        iter = iter,
                        refresh = refresh)
    )
    expect_geostan(fit)
})


test_that("CAR accepts covariate ME with WX, mixed ME-non-ME", {
    data(georgia)
    A <- shape2mat(georgia)
    cars <- prep_car_data(A)
    ME <- prep_me_data(se = data.frame(insurance = georgia$insurance.se),
                 bounds = c(0, 100),
                 car_parts = cars)
    SW(
        fit <- stan_car(log(rate.male) ~ insurance + ICE,
                        slx = ~ insurance + ICE,
                        data = georgia,
                        ME = ME,
                        car_parts = cars,
                        chains = 1,
                        iter = iter,
                        refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("DCAR example runs", {
    A <- shape2mat(georgia, "B")
    D <- sf::st_distance(sf::st_centroid(georgia))
    A <- D * A
    cp <- prep_car_data(A, "DCAR", k = 1)
    fit <- stan_car(log(rate.male) ~ college,
                    data = georgia,
                    car = cp,
                    iter = iter,
                    chains = 1)
     expect_geostan(fit)
})

test_that("CAR with censored y", {
    data(georgia)
    A <- shape2mat(georgia)
    cars <- prep_car_data(A)    
    SW(
        fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)) + ICE + college,
                        censor_point = 9,
                        data = georgia,
                        chains = 1,
                        family = poisson(),
                        car_parts = cars,
                        iter = iter,
                        refresh = refresh
                        )
    )
    expect_geostan(fit)
})

test_that("Slim SAR works", {
    row = 20
    col = 26
    N <- row * col
    cdl <- prep_car_data2(row = row, col = col)
    x <- rnorm(n = N)
    y <- .75 *x + rnorm(n = N, sd = .5)
    df <- data.frame(y=y, x=x)
    fit <- stan_car(y ~ x,
                    slx = ~ x,
                    data = df,
                    car_parts = cdl,
                    chains = 1,
                    iter = iter,
                    refresh = refresh,
                    slim = TRUE)
    expect_geostan(fit)
})
