
iter=10
refresh=0
source("helpers.R")

context("stan_car")
test_that("Poisson offset model works, car", {
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    ME <- list(offset = rep(10, n))
    SW(
        fit <- stan_car(sents ~ offset(expected_sents),
                    data = sentencing,
                    ME = ME,
                    C = C,
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)    
})

test_that("CAR accepts covariate ME, mixed (un-) bounded", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.75, n),
                          historic_gop = rep(3, n)),
               bounded = c(1, 0))
    SW(
        fit <- stan_car(cbind(trump_2016, total_2016 - trump_2016) ~ unemployment + historic_gop,
                        data = ohio,
                        C = C,
                        family = binomial(),
                        ME = ME,
                        chains = 1,
                        iter = iter,
                        refresh = refresh)
    )
    expect_geostan(fit)
})



test_that("CAR accepts covariate ME with WX, mixed ME-non-ME", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.75, n),
                               historic_gop = rep(3, n),
                               college_educated = rep(3, n)),
               bounded = c(1, 0, 1))
    SW(
        fit <- stan_car(cbind(trump_2016, total_2016 - trump_2016) ~ log(population) + college_educated + unemployment + historic_gop,
                    slx = ~ college_educated + unemployment + log(population),
                    data = ohio,
                    family = binomial(),
                    C = C,
                    ME = ME,
                    chains = 1,
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})
