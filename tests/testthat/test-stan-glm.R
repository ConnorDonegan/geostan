iter=20
refresh = 0
source("helpers.R")

context("stan_glm")
test_that("GLM works", {
    data(sentencing)
    SW(
        fit <- stan_glm(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    chains = 1,
                    family = poisson(),
                    init_r = 0.2,
                    iter = iter,
                    refresh = refresh)
       )
    expect_geostan(fit)
    SW(
        fit <- stan_glm(log(sents/expected_sents) ~ 1,
                    data = sentencing,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
       )
    expect_geostan(fit)    
})

test_that("GLM works, centerx = TRUE", {
    SW(
        fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + college + black,
                    data = georgia,
                    chains = 1,
                    centerx = TRUE,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
       )
    expect_geostan(fit)
})


test_that("GLM works with covariate ME", {
    data(georgia)
    n <- nrow(georgia)
    ME <- prep_me_data(se = data.frame(ICE = georgia$ICE.se))
    SW(
        fit <- stan_glm(log(rate.male) ~ ICE,
                    ME = ME,                        
                    data = georgia,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("GLM works with covariate ME: spatial data model", {
    data(georgia)
    A <- shape2mat(georgia, "B")
    ME <- prep_me_data(
        se = data.frame(ICE = georgia$ICE.se),
        car_parts = prep_car_data(A)
        )
    SW(
        fit <- stan_glm(log(rate.male) ~ ICE,
                    ME = ME,                        
                    data = georgia,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("GLM accepts covariate ME, multiple proportions", {
    data(georgia)
    A <- shape2mat(georgia, "B")
    ME <- prep_me_data(
        se = data.frame(
            insurance = georgia$insurance.se,
            college = georgia$college.se
            ),
        bounds = c(0, 100)
        )
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance + college,
                    ME = ME,                        
                    data = georgia,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("GLM accepts covariate ME, with logit transform", {
    data(georgia)
    georgia$insurance <- georgia$insurance / 100
    georgia$insurance.se <- georgia$insurance.se / 100
    ME <- prep_me_data(
        se = data.frame(
            insurance = georgia$insurance.se
        ),
        logit = TRUE,
        bounds = c(0, 1)
        )
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance,
                    ME = ME,                        
                    data = georgia,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance + college,
                    ME = ME,                        
                    data = georgia,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})


test_that("Set priors for GLM", {
    data(georgia)
    SW(
        fit <- stan_glm(log(rate.male) ~ insurance + college,
                        data = georgia,
                        chains = 1,
                        family = student_t(),
                        prior = list(
                            nu = gamma(3, 0.1),
                            beta = normal(location = c(0,0),
                                                   scale = c(10,10)),
                                     intercept = normal(0, 10),
                                     sigma = student_t(10, 0, 5)
                                     ),
                    iter = iter,
                    refresh = refresh,
                    init_r = 0.1
                    )
    )
    expect_geostan(fit)
})



test_that("GLM with censored y", {
    data(georgia)
    SW(
        fit <- stan_glm(deaths.female ~ offset(log(pop.at.risk.female)) + ICE + college,
                        censor_point = 9,
                        data = georgia,
                        chains = 1,
                        family = poisson(),
                        prior = list(
                            beta = normal(location = c(0,0),
                                                   scale = c(5,5)),
                                     intercept = normal(-4, 5)
                                     ),
                    iter = iter,
                    refresh = refresh,
                    init_r = 0.1
                    )
    )
    expect_geostan(fit)
})


test_that("Slim GLM works", {
    N <- 20 * 26
    x <- rnorm(n = N)
    y <- .75 *x + rnorm(n = N, sd = .5)
    df <- data.frame(y=y, x=x)
    fit <- stan_glm(y ~ x,
                    data = df,
                    chains = 1,
                    iter = iter,
                    refresh = refresh,
                    slim = TRUE)
    expect_geostan(fit)
})
