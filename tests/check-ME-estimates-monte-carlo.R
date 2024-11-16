
##
## A Monte Carlo study to test measurement error models
##

library(geostan)

## no. iterations
M = 15

## a regular grid
ncol = 20
nrow = 20
sars <- prep_sar_data2(ncol, nrow, quiet = TRUE)
cars <- prep_car_data2(ncol, nrow, quiet = TRUE)
W <- sars$W
N <- nrow(W)

# iterate:
b = 0.5
g = -0.5
rho.x = 0 ## checking ME, not checking spatial models
rho.y = 0.6
sigma.y = 0.3
sigma.me = 0.3
pars <- c(const = 0,
          gamma = g,
          beta = b,
          rho = rho.y,
          scale = sigma.y)

res <- sapply(1:M, FUN = function(i) {

    x <- sim_sar(w=W, rho=rho.x)
    Wx <- (W %*% x)[,1]    
    mu <- b * x + g * Wx    
    y <- sim_sar(w=W, rho=rho.y, mu = mu, sigma = sigma.y, type = "SEM")
    x = x + rnorm(N, sd = sigma.me)        
    dat <- data.frame(y, x)    

    ME <- prep_me_data(se = data.frame(x = rep(sigma.me, N)))
    if (rho.x > 0) ME <- prep_me_data(se = data.frame(x = rep(sigma.me, N)), car_parts = cars)

    fit_sem_me <- geostan::stan_sar(y ~ x,
                                 data = dat,
                                 type = "SDEM",
                                 sar_parts = sars,
                                 ME = ME,
                                 chains = 1,
                                 iter = 600,
                                 quiet = TRUE,
                                 slim = TRUE
                                 ) |>
        suppressWarnings()

    fit_sem <- geostan::stan_sar(y ~ x,
                                  data = dat,
                                  type = "SDEM",
                                  sar_parts = sars,
                                  ## ME = ME, ##
                                  chains = 1,
                                  iter = 600,
                                  quiet = TRUE,
                                  slim = TRUE
                                  ) |>
        suppressWarnings()
    ## fit_glm0 <- geostan::stan_glm(y ~ x,
        ##                          slx = ~ x,
        ##                          data = dat,
        ##                          ## ME = ME, ##
        ##                          C = cars$C,
        ##                          chains = 1,
        ##                          iter = 600,
        ##                          quiet = TRUE,
        ##                          slim = TRUE
        ##                          ) |>
        ## suppressWarnings()
    ## fit_car <- geostan::stan_car(y ~ x,
    ##                              slx = ~ x,
    ##                              data = dat,
    ##                              car_parts = cars,
    ##                              ME = ME,
    ##                              chains = 1,
    ##                              iter = 900,
    ##                              quiet = TRUE,
    ##                              slim = TRUE
    ##                              ) |>
    ##     suppressWarnings()

    ## fit_glm <- geostan::stan_glm(y ~ x,
    ##                              slx = ~ x,
    ##                              data = dat,
    ##                              ME = ME,
    ##                              C = cars$C,
    ##                              chains = 1,
    ##                              iter = 900,
    ##                              quiet = TRUE,
    ##                              slim = TRUE
    ##                              ) |>
    ##     suppressWarnings()
    
    x <- c(
        # GLM0 = fit_glm0$summary[c('intercept', 'w.x', 'x', 'rho', 'sigma'), 'mean'],
       # GLM = fit_glm$summary[c('intercept', 'w.x', 'x', 'rho', 'sigma'), 'mean'],        
        SEM_ME = fit_sem_me$summary[c('intercept', 'w.x', 'x', 'sar_rho', 'sar_scale'), 'mean'],
        SEM = fit_sem$summary[c('intercept', 'w.x', 'x', 'sar_rho', 'sar_scale'), 'mean'] #,        
       # CAR = fit_car$summary[c('intercept', 'w.x', 'x', 'car_rho', 'car_scale'), 'mean']
        )
    
    return (x)
    
    })


RMSE <- function(est, true) sqrt(mean(est - true)^2)

g_res <- res |>
    subset(row.names(res) %in% c('SEM_ME2', 'SEM2'))
   ## subset(row.names(res) %in% c('GLM02', 'GLM2', 'SEM2', 'CAR2'))
b_res <- res |>
    subset(row.names(res) %in% c('SEM_ME3', 'SEM3'))    
     ##  subset(row.names(res) %in% c('GLM03', 'GLM2', 'SEM3', 'CAR3'))

cat("\n**\nRMSE of ME models: \n**\n")

cat("\n**\nGamma\n**\n")
apply(g_res, 1, RMSE, g) |>
    print(digits = 2)

cat("\n**\nBeta\n**\n")
apply(b_res, 1, RMSE, b) |>
    print(digits = 2)


cat("\n**\nME models\nSEM Estimates (rounded) should be close to DGP parameters \n",M, "iterations\n**\n")

est <- apply(res, 1, mean) 

est <- est[c('SEM2', 'SEM3', 'SEM_ME2', 'SEM_ME3')]

out <- data.frame(
   # Mod = attributes(est)$names,
    Par = rep(c('Gamma', 'Beta'), 2),
    DGP = c(g, b, g, b),
    Est = est
) |>
    transform(
        Error = round(Est - DGP, 2)
    ) 

out |>
    print()

write.csv(out, "check-ME-estimates-monte-carlo-output.csv")
