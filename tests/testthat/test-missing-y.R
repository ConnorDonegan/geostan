iter=10
refresh = 0
source("helpers.R")

context("missing y")

test_that("Poisson missing y", {
    data(georgia)
    N = nrow(georgia)
  georgia$deaths.female[sample.int(N, size = 12)] <- NA
  georgia$deaths.female[sample.int(N, size = 3)] <- 0
  SW(
    fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)),
                    data = georgia,
                    car_parts = prep_car_data(shape2mat(georgia, "B")),
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
  )
  expect_geostan(fit)    
})

test_that("Binomial missing y", {
    data(georgia)
    A = shape2mat(georgia, "B")
    N = nrow(A)
    georgia$deaths.female[sample.int(N, size = 25)] <- NA    
    georgia$y <- round(georgia$deaths.female / 10)
    georgia$y[sample.int(N, 5)] <- 0
    georgia$f <- round(4 * georgia$deaths.female / 10)
    georgia$f[which(!is.na(georgia$y))[1]] <- NA
  SW(
    fit <- stan_glm(cbind(y, f) ~ 1,
                    data = georgia,
                  #  car_parts = prep_car_data(shape2mat(georgia, "B")),
                    chains = 1,
                    family = binomial(),
                    iter = iter,
                    refresh = refresh)
  )
  expect_geostan(fit)
    SW(
        
    fit <- stan_icar(cbind(y, f) ~ 1,
                    data = georgia,
                    C = A,
                    chains = 1,
                    family = binomial(),
                    iter = iter,
                    refresh = refresh)
    
  )

})


test_that("ESF missing y", {
  data(georgia)
  georgia$deaths.female[1:10] <- NA
  georgia$y <- georgia$deaths.female
  georgia$f <- round(4 * georgia$deaths.female)
  SW(
  fit <- stan_esf(cbind(y, f) ~ log(income),
                    data = georgia,
                  C = shape2mat(georgia, "B"),
                    chains = 1,
                    family = binomial(),
                    iter = iter,
                    refresh = refresh)
  )
  expect_geostan(fit)    
})



