iter=30
refresh = 0
source("helpers.R")

context("stan_esf")

test_that("make_EV returns a data.frame", {
    A <- shape2mat(georgia, "B")
    EV <- make_EV(A)
    expect_is(EV, "data.frame")
})


test_that('data types correct', {
    A <- shape2mat(georgia, "B")
    expect_is(A, 'Matrix')
    EV <- make_EV(A)
    expect_is(EV, 'data.frame')
    EV.l <- make_EV(A, values = TRUE)
    expect_is(EV.l, 'list')
    L <- lisa(EV[,10], A)
    expect_is(L, 'data.frame')
    L.df <- lisa(EV[,10], A, type = FALSE)
    expect_is(L.df, 'numeric')
})

test_that("ESF model works", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    SW(
        fit <- stan_esf(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
       )
    expect_geostan(fit)
})

test_that("ESF model works by providing EV", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing, "B")
    EV <- make_EV(C)
    SW(
        fit <- stan_esf(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    EV = EV,
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
                fit <- stan_esf(deaths.male ~ offset(log(pop.at.risk.male)) + college + log_income,
                                data = georgia,
                                re = ~ GEOID,
                                 ME = ME,
                                 C = C,                                
                                 prior = list(
                                     intercept = normal(0, 5),
                                     beta = normal(c(0,0), c(4,4))
                                 ),
                                 chains = 1,
                                 family = poisson(),
                                 iter = iter,
                                 refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("ESF with censored y", {
    data(georgia)
    C <- shape2mat(georgia)
    SW(
        fit <- stan_esf(deaths.female ~ offset(log(pop.at.risk.female)) + insurance,
                        re = ~ GEOID,
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
