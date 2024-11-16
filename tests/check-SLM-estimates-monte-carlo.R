

##
## verify estimates of key model parameters for SDLM, including impacts
##

#devtools::load_all("~/dev/geostan")
#library(spatialreg)
library(geostan)

# no. iterations [s > 90 required for impact est. to converge - they are high variance; others converge by 30]
S <- 50

# use measurement error in x
has_me <- FALSE

# regular lattice
parts <- prep_sar_data2(row = 10, col = 30, quiet = TRUE)
cp = prep_car_data2(10, 30, quiet=T)
W <- parts$W

# for spatialreg
#trm <- spatialreg::trW(W, type="mult")
#lw <- spdep::mat2listw(W, style = 'W')

# model parameters
N <- nrow(W)
B <- -0.5
G <- 0.1
R <- 0.7
Sig = 0.1
sigma_me <- 0.1
pars <- c(const = 0, beta = B, gamma = G, rho = R, sigma = Sig)
pars <- c(pars, geostan::spill(beta = B, gamma = G, rho = R, W = W, approx=FALSE))

# results for S iterations
## res <- res2 <- matrix(NA, nrow = S, ncol = 8)
res <- matrix(NA, nrow = S, ncol = 8)

for (s in 1:S) {
    x <- sim_sar(w=W, rho=R)
    if (has_me) x = x + rnorm(N, sd = sigma_me)
    Wx <- (W %*% x)[,1]
    mu <- B * x + G * Wx    
    y <- sim_sar(w=W, rho=R, mu = mu, sigma = Sig, type = "SLM")
    dat <- data.frame(y, x)    

    ME <- prep_me_data(se = data.frame(x = rep(sigma_me, N)))
    
    if (has_me) {
        fit <- stan_sar(y ~ x,
                        data = dat,
                        sar = parts,
                        ME = ME,
                        type = "SDLM",
                        iter = 400,
                        chains = 1,
                        slim = TRUE,
                        quiet = TRUE) |>
            suppressWarnings()
    } else {
        fit <- stan_sar(y ~ x,
                        data = dat,
                        sar = parts,
                        type = "SDLM",
                        iter = 400,
                        chains = 1,
                        slim = TRUE,
                        quiet = TRUE) |>
            suppressWarnings()
    }
    
    res[s, ] <- c(
        fit$summary[c('intercept', 'x', 'w.x', 'sar_rho', 'sar_scale'), 'mean'],
        geostan::impacts(fit, approx = FALSE)$summary$x[,'mean']
    )

    ## fit2 <- spBreg_lag(y ~ x, data = dat, listw=lw, Durbin = TRUE,
    ##                    control =  list(ndraw = 1000L))
    ## x <- spatialreg::impacts(fit2, tr = trm)$sres
    ## res2[s, ] <- c(
    ##     apply(fit2, 2, mean),
    ##     mean(x$direct),
    ##     mean(x$indirect),
    ##     mean(x$total)
    ## )
}

# spatialreg reports variance parameter, geostan reports scale param.
## res2[, 5] <- sqrt(res2[,5])

cat("\n**\nSDLM Monte Carlo results\n", S, "iterations\nAverage values\n**\n")

est <- apply(res, 2, mean)

dat <- cbind(
    DGP = pars,
    geostan =  est,
    Error = est - pars
      #spatialreg = apply(res2, 2, mean)
) |>
    round(2)

dat |>
    print()

write.csv(dat, "check-SLM-estimates-monte-carlo-output.csv")


## # viz comparison of all estimates
## par(mfrow = c(2, 4), mar = c(4, 4, 1, 1))
## for (j in 1:8) {
##     plot(res[, j], res2[, j],
##          main = attributes(pars)$names[j],
##          xlab = NA,
##          ylab = NA,
##          bty='L', pch=22);
##     abline(0,1, lty = 3, lwd = .5)
##     abline(v = pars[j], col = rgb(0.1, 0.3, 0.5, 0.5))
##     abline(h = pars[j], col = rgb(0.1, 0.3, 0.5, 0.5))
##     mtext("geostan", side = 1, line = 2)
##     mtext("spatialreg", side = 2, line = 2)    
## }


