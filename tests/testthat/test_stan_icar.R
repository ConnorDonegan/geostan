
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

test_that("IAR accepts covariate ME, multiple bounded x vars", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 1))
    SW(
        fit <- stan_icar(cbind(trump_2016, total_2016 - trump_2016) ~ unemployment + historic_gop,
                        data = ohio,
                        C = C,
                        family = binomial(),
                        ME = ME,
                        chains = 1,
                        iter = iter,
                        silent = silent)
       )
    expect_geostan(fit)
})

test_that("IAR accepts covariate ME, mixed (un-) bounded", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 0))
    SW(
        fit <- stan_icar(cbind(trump_2016, total_2016 - trump_2016) ~ unemployment + historic_gop,
                        data = ohio,
                        C = C,
                        family = binomial(),
                        ME = ME,
                        chains = 1,
                        iter = iter,
                        silent = silent)
    )
    expect_geostan(fit)
})


test_that("IAR accepts covariate ME with WX, mixed ME-non-ME", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.75, n),
                               historic_gop = rep(3, n),
                               college_educated = rep(3, n)),
               bounded = c(1, 0, 1))
    SW(
        fit <- stan_icar(cbind(trump_2016, total_2016 - trump_2016) ~ log(population) + college_educated + unemployment + historic_gop,
                    slx = ~ college_educated + unemployment + log(population),
                    data = ohio,
                    family = binomial(),
                    C = C,
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    silent = silent)
    )
    expect_geostan(fit)
})


