
iter=10
refresh = 0
source("helpers.R")

context("stan_icar")
test_that("Poisson model works, icar: icar", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    SW(
        fit <- stan_icar(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    type = "icar",
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)    
})


test_that("Poisson offset model works, icar: bym", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    SW(
        fit <- stan_icar(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    type = "bym",
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)    
})

test_that("Poisson model works, icar: bym2", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    SW(
        fit <- stan_icar(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    type = "bym2",
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)    
})


test_that("Poisson model works with length-2 vector beta prior", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    SW(
        fit <- stan_icar(sents ~ offset(log(expected_sents)) + plantation_belt,
                    data = sentencing,
                    type = "bym2",
                    prior = list(
                        intercept = normal(0, 5),
                        beta = normal(0, 2)
                        ),
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)    
})

test_that("ICAR accepts covariate ME with logit transform", {
    data(georgia)
    C <- shape2mat(georgia)    
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
                fit <- stan_icar(deaths.male ~ offset(log(pop.at.risk.male)) + college + log_income,
                                 data = georgia,
                                 ME = ME,
                                 type = "bym2",
                                 prior = list(
                                     intercept = normal(0, 5),
                                     beta = normal(c(0,0), c(4,4))
                                 ),
                                 C = C,
                                 chains = 1,
                                 family = poisson(),
                                 iter = iter,
                                 refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("BYM with censored y", {
    data(georgia)
    C <- shape2mat(georgia)
    SW(
        fit <- stan_icar(deaths.female ~ offset(log(pop.at.risk.female)),
                         re = ~ GEOID,
                         type = "bym",
                        censor_point = 9,
                        data = georgia,
                        chains = 1,
                        family = poisson(),
                        C = C,
                        iter = iter,
                        refresh = refresh
                        )
    )
    expect_geostan(fit)
})
