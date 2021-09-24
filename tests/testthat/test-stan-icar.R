
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
