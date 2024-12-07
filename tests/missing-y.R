
##
## models with missing outcome data
##

#devtools::load_all("~/dev/geostan")
library(geostan)

iter = 25

# CAR Poisson model
data(georgia)
C <- shape2mat(georgia, "B", quiet = TRUE)
cars <- prep_car_data(C, quiet = TRUE)
N = nrow(georgia)
georgia$deaths.female[sample.int(N, size = 12)] <- NA
georgia$deaths.female[sample.int(N, size = 3)] <- 0

fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)),
                data = georgia,
                car_parts = cars,
                chains = 1,
                family = poisson(),
                iter = iter,
                quiet = TRUE) |>
    suppressWarnings()


# Binomial model with missing y
data(georgia)
A = shape2mat(georgia, "B")
N = nrow(A)
georgia$deaths.female[sample.int(N, size = 25)] <- NA    
georgia$y <- round(georgia$deaths.female / 10)
georgia$y[sample.int(N, 5)] <- 0
georgia$f <- round(4 * georgia$deaths.female / 10)


# glm
fit <- stan_glm(cbind(y, f) ~ 1,
                data = georgia,
                chains = 1,
                family = binomial(),
                iter = iter,
                quiet = TRUE) |>
    suppressWarnings()

# icar
fit <- stan_icar(cbind(deaths.female, pop.at.risk.female) ~ 1,
                 data = georgia,
                 type = 'bym',
                 C = A,
                 chains = 1,
                 family = binomial(),
                 iter = iter,
                 quiet = TRUE) |>
    suppressWarnings()

# esf
data(georgia)
georgia$deaths.female[1:10] <- NA
georgia$y <- georgia$deaths.female
georgia$f <- round(4 * georgia$deaths.female)

fit <- stan_esf(cbind(y, f) ~ log(income),
                data = georgia,
                C = shape2mat(georgia, "B"),
                chains = 1,
                family = binomial(),
                iter = iter,
                quiet = TRUE) |>
    suppressWarnings()

