
iter=10
silent = TRUE
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
                    silent = silent)
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
                    silent = silent)
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
                    silent = silent)
    )
    expect_geostan(fit)    
})

