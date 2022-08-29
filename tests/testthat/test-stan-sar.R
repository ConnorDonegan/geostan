iter=10
refresh = 0
source("helpers.R")

context("stan_sar")

test_that("Poisson SAR model works", {
    data(sentencing)
    SW(
        fit <- stan_sar(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    C = shape2mat(sentencing, "W"),
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
  )
    expect_geostan(fit)    
})

test_that("SAR accepts covariate ME with logit transform", {
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


