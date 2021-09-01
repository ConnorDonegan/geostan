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
    expect_is(A, 'matrix')
    EV <- make_EV(A)
    expect_is(EV, 'data.frame')
    EV.l <- make_EV(A, values = TRUE)
    expect_is(EV.l, 'list')
    L <- lisa(EV[,10], A)
    expect_is(L, 'numeric')
    L.df <- lisa(EV[,10], A, type = TRUE)
    expect_is(L.df, 'data.frame')
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




