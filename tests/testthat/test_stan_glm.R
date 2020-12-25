iter=20
silent = TRUE
source("helpers.R")

context("stan_glm")
test_that("Poisson offset model works", {
    data(sentencing)
    n <- nrow(sentencing)
    ME <- list(offset = rep(10, n))
    SW(
        fit <- stan_glm(sents ~ offset(log(expected_sents)),
                    data = sentencing,
                    ME = ME,
                    chains = 1,
                    family = poisson(),
                    init_r = 0.2,
                    iter = iter,
                    silent = silent)
       )
    expect_geostan(fit)
})

test_that("GLM works with covariate ME", {
    data(ohio)
    n <- nrow(ohio)
    ME <- list(se = data.frame(unemployment = rep(0.75, n)))
    SW(
        fit <- stan_glm(gop_growth ~ unemployment + historic_gop,
                    data = ohio,
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})

test_that("GLM works with covariate ME: spatial data model", {
    data(ohio)
    C <- shape2mat(ohio, "B")
    n <- nrow(ohio)
    ME <- list(se = data.frame(unemployment = rep(0.75, n)),
               spatial = TRUE
               )
    SW(
        fit <- stan_glm(gop_growth ~ unemployment + historic_gop,
                    data = ohio,
                    ME = ME,
                    C = C,
                    chains = 1,
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})

test_that("GLM accepts covariate ME, multiple x proportions", {
    data(ohio)
    n <- nrow(ohio)
    ME <- list(se = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 1))
    SW(
        fit <- stan_glm(gop_growth ~ unemployment + historic_gop,
                    data = ohio,
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    silent = silent)
       )
    expect_geostan(fit)
})

test_that("GLM accepts covariate ME, mixed (un-) bounded with user-defined bounds", {
    data(ohio)
    n <- nrow(ohio)
    ohio$unemployment <- ohio$unemployment / 100
    ME <- list(se = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 0),
               bounds = c(0, 1))
    SW(
        fit <- stan_glm(gop_growth ~ unemployment + historic_gop,
                    data = ohio,
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})

test_that("GLM accepts covariate ME with WX, mixed ME-non-ME", {
    data(ohio)
    n <- nrow(ohio)
    ME <- list(se = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 0))
    SW(
        fit <- stan_glm(gop_growth ~ log(population) + college_educated + unemployment + historic_gop,
                    slx = ~ college_educated + unemployment,
                    data = ohio,
                    C = shape2mat(ohio),
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})

test_that("Binomial GLM accepts covariate ME with WX, mixed ME-non-ME, no bounded term", {
    data(ohio)
    n <- nrow(ohio)
    ME <- list(se = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)))
    SW(
        fit <- stan_glm(cbind(trump_2016, total_2016 - trump_2016) ~ log(population) + college_educated + unemployment + historic_gop,
                    slx = ~ college_educated + unemployment,
                    data = ohio,
                    C = shape2mat(ohio),
                    ME = ME,
                    chains = 1,
                    family = binomial(),
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})




