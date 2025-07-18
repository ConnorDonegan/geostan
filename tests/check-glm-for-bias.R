
# this should not run on CRAN or in every test instance.

RUN_ME = FALSE
if (RUN_ME) {
    
###########################################
## this line overwrites previous results ##
out_file <- "check-glm-for-bias"
writeLines(paste0("# Results from tests/check-glm-for-bias.R \n# Run date: ", Sys.Date()),
           out_file)
###########################################

library(geostan)
library(parallel)
nsim <- 2
pars <- c("alpha", "beta1", "beta2")
W <- shape2mat(georgia, "W")
A <- shape2mat(georgia, "B")
P <- georgia$population
N <- nrow(georgia)
CP <- prep_car_data(A)
EV <- make_EV(A)

    
sim_params <- function() {    
        rho_x <- rbeta(n = 1, 3, 3)
        rho_phi <- rbeta(n = 1, 3, 3)
        beta <- runif(n = 2, -2, 2)
        alpha <- truncnorm::rtruncnorm(n = 1, mean = -5.5, sd = 1, a = -Inf, b = 0)
        sigma_phi <- 0.3
        l <- list(rho_x = rho_x,
             rho_phi = rho_phi,
             beta = beta,
             alpha = alpha,
             sigma_phi = sigma_phi
             )
        return (l)
}

sim_xy <- function(rho_x, alpha, beta, rho_phi, sigma_phi) {
        x <- sim_sar(m = 2, w = W, rho= rho_x, sigma = 0.5)
        x = scale(t(x), center = TRUE, scale = FALSE)
        phi <- sim_sar(m = 1, w = W, mu = as.numeric(alpha + x %*% beta), sigma = sigma_phi, rho = rho_phi)
        y <- rpois(n = N, lambda = exp(phi) * P)
        df <- data.frame(y = y, x1 = x[,1], x2 = x[,2], P = P, id = as.character(1:length(y)))
        return(df)
}

message("testing that CAR: Log-linear spatial Poisson SLX models are unbiased")
# temp file for results
tmp_car <- tempfile(fileext = "txt")
    # simulate spatial x, spatial y; fit model y ~ wx + x
fit_fn_car <- function(i) {
    plist <- sim_params()
    df <- sim_xy(plist$rho_x, plist$alpha, plist$beta, plist$rho_phi, plist$sigma_phi)
                                        # CAR
    fit_car <- stan_car(y ~ offset(log(P)) + x1 + x2,
                       # slx = ~ x1 + x2,
                        data = df,
                        family = poisson(),
                        car_parts = CP,
                        chains = 1,
                        iter = 1200
                        )$summary
    res_car <- c(fit_car[c('intercept', 'x1', 'x2'), 'mean'] - c(plist$alpha, plist$beta)) 
    write.table(matrix(res_car,nrow=1), file = tmp_car, append = TRUE, row.names = FALSE, col.names = FALSE)
    cat(i, " (CAR) ")
}

# run nsim times
quietly=mclapply(1:nsim, fit_fn_car, mc.cores = 5)
## mean errors should all be zero
res_car <- read.table(tmp_car)
bias_car <- apply(res_car, 2, mean)
std_error <- apply(res_car, 2, sd) / sqrt(nsim)
#for (i in 1:length(bias_car)) testthat::expect_equal(as.numeric(bias_car[i]), 0, tol = 2 * std_error)
cat("# Context: CAR log-linear Poisson SLX models are unbiased", file = out_file, append = TRUE)
cat(paste0("\n# ", pars, " ", round(bias_car, 3), " Expected result: 0 \\pm ", round(2 * std_error, 3)),
    file = out_file,
    append = TRUE)

message("Testing that BYM: Log-linear spatial Poisson SLX models are unbiased")
    # temp file for results
tmp_icar <- tempfile(fileext = "txt")
    # simulate spatial x, spatial y; fit model y ~ wx + x
fit_fn_icar <- function(i) {
    plist <- sim_params()
    df <- sim_xy(plist$rho_x, plist$alpha, plist$beta, plist$rho_phi, plist$sigma_phi)
        # ESF
    fit_icar <- stan_icar(y ~ offset(log(P)) + x1 + x2,
                          type = "bym",
                          data = df,
                          family = poisson(),
                          C = A,
                          chains = 1,
                          iter = 1500
                          )$summary
    res_icar <- c(fit_icar[c('intercept', 'x1', 'x2'), 'mean'] - c(plist$alpha, plist$beta)) 
    write.table(matrix(res_icar,nrow=1), file = tmp_icar, append = TRUE, row.names = FALSE, col.names = FALSE)
    cat(i, " (ICAR)")
}       
quietly=mclapply(1:nsim, fit_fn_icar, mc.cores = 5)
## mean errors should all be zero
res_icar <- read.table(tmp_icar)
bias_icar <- apply(res_icar, 2, mean)
std_error <- apply(res_icar, 2, sd) / sqrt(nsim)    
cat("\n# Context: BYM log-linear Poisson SLX models are unbiased", file = out_file, append = TRUE)
cat(paste0("\n# ", pars, " ", round(bias_icar, 3), " Expected result: 0 \\pm ", round(2 * std_error, 3)),
    file = out_file,
    append = TRUE)

 
message("Testing that ESF: Log-linear spatial Poisson SLX models are unbiased")
    # temp file for results
tmp_esf <- tempfile(fileext = "txt")
                                        # simulate spatial x, spatial y; fit model y ~ wx + x    
fit_fn_esf <- function(i) {
    plist <- sim_params()
    df <- sim_xy(plist$rho_x, plist$alpha, plist$beta, plist$rho_phi, plist$sigma_phi)
        # ESF
    fit_esf <- stan_esf(y ~ offset(log(P)) + x1 + x2,
                        re = ~ id,
                        data = df,
                        family = poisson(),
                        EV = EV,
                        C = A,
                        chains = 1,
                        iter = 1500
                        )$summary
    res_esf <- c(fit_esf[c('intercept', 'x1', 'x2'), 'mean'] - c(plist$alpha, plist$beta)) 
    write.table(matrix(res_esf,nrow=1), file = tmp_esf, append = TRUE, row.names = FALSE, col.names = FALSE)
    cat(i, " (ESF)")
}

quietly=mclapply(1:nsim, fit_fn_esf, mc.cores = 5)
## mean errors should all be zero
res_esf <- read.table(tmp_esf)
bias_esf <- apply(res_esf, 2, mean)
std_error <- apply(res_esf, 2, sd) / sqrt(nsim)
cat("\n# Context: ESF log-linear Poisson SLX models are unbiased", file = out_file, append = TRUE)
cat(paste0("\n# ", pars, " ", round(bias_esf, 3), " Expected result: 0 \\pm ", round(2 * std_error, 3)),
    file = out_file,
    append = TRUE)
}
