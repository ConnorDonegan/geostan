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

test_that("GLM works with covariate ME", {
    data(georgia)
    n <- nrow(georgia)
    ME <- list(se = data.frame(ICE = georgia$ICE.se))
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
    ME <- list(
        se = data.frame(ICE = georgia$ICE.se),
        spatial = TRUE,
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
    ME <- list(
        se = data.frame(
            insurance = georgia$insurance.se,
            college = georgia$college.se
            ),
        spatial = FALSE,
        bounded = c(1,1),
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


test_that("Set priors for GLM", {
    data(georgia)
    SW(
        fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + insurance + college,
                        data = georgia,
                        chains = 1,
                        family = poisson(),
                        prior = list(beta = data.frame(location = c(0,0),
                                                       scale = c(10,10)),
                                     intercept = c(0, 10)
                                     ),
                    iter = iter,
                    refresh = refresh,
                    init_r = 0.1
                    )
    )
    expect_geostan(fit)
})

test_that("Set length-2 vector prior for 1 beta in GLM", {
    data(georgia)
    SW(
        fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + insurance,
                        data = georgia,
                        chains = 1,
                        family = poisson(),
                        prior = list(beta = c(0, 10),
                                     intercept = c(0, 10)
                                     ),
                    iter = iter,
                    refresh = refresh,
                    init_r = 0.1
                    )
    )
    expect_geostan(fit)
})

