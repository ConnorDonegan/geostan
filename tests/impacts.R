




library(geostan)
library(spatialreg)


stop("THERE IS A PROBLEM WITH SLX.")



parts <- prep_sar_data2(row = 10, col = 30, quiet = TRUE)
W <- parts$W

trm <- spatialreg::trW(W, type="mult")
lw <- spdep::mat2listw(W, style = 'W')

N <- nrow(W)
S <- 50
B <- -0.5
G <- 0.1
R <- 0.6
sigma = 0.5
pars <- c(0, B, G, R, 1)
res <- res2 <- matrix(NA, nrow = S, ncol = 4)#5, 8
for (s in 1:S) {
    x <- sim_sar(w=W, rho=R)
    Wx <- (W %*% x)[,1]
    mu <- B * x + G * Wx    
    y <- sim_sar(w=W, rho=R, mu = mu, sigma = sigma)
    dat <- data.frame(y, x, Wx)
    ME <- prep_me_data(se = data.frame(x = rep(0.01, N)))

    mat <- cbind(Wx, x)
    qm <- qr(mat)

    XB <- mat %*% c(G, B)
    XB2 <- qr.Q(qm) %*% qr.R(qm) %*% c(G, B)    

    qr.R(qm)[,1] %*% matrix(G)
    qr.R(qm) %*% c(G, B)
    
    fit <- stan_glm(y ~  x,
                    slx = ~ x,
                    ME = ME,
                    data = dat,
                    C = W,
                    centerx = TRUE,
                    iter = 600,
                    chains = 1,
                    slim = TRUE,
                    quiet = TRUE) |>
        suppressWarnings()
    res[s, ] <- c(
        fit$summary[c('intercept', 'x', 'w.x', 'sigma'), 'mean']#,
        #geostan::impacts(fit)$summary$x[,'mean']
    )
}

apply(res, 2, mean) |> round(3)


    ## fit <- stan_sar(y ~ x,
    ##                 data = dat,
    ##                 sar = parts,
    ##                 type = "SDEM",
    ##                 iter = 600,
    ##                 chains = 1,
    ##                 slim = TRUE,
    ##                 quiet = TRUE) |>
    ##     suppressWarnings()
    
    res[s, ] <- c(
        fit$summary[c('intercept', 'x', 'w.x', 'sar_rho', 'sar_scale'), 'mean']#,
        #geostan::impacts(fit)$summary$x[,'mean']
    )

    fit2 <- spBreg_err(y ~ x, data = dat, listw=lw, Durbin = TRUE,
                       control =  list(ndraw = 1e3L))
    #x <- spatialreg::impacts(fit2, tr = trm)$sres
    res2[s, ] <- c(
        apply(fit2, 2, mean)#,
        #mean(x$direct),
        #mean(x$indirect),
        #mean(x$total)
    )
}


apply(res, 2, mean) |> round(1)
apply(res2[, 1:5], 2, mean) |> round(1)

par(mfrow = c(2, 4), mar = rep(1,4))
for (j in 1:8) {
    plot(res[, j], res2[, j], bty='L', pch=22);
    abline(0,1)
}






colMeans(res)

spill(0.5, 0.25, 0.6, W)

og = par(mfrow = c(2, 3),
         mar = c(3, 3, 1, 1))
hist(res[,1], main = 'Direct'); abline(v=0.611, col = 'green')
hist(res[,2], main = 'Indirect'); abline(v=1.26, col='green')
hist(res[,3], main = 'Total'); abline(v=1.87, col='green')
hist(res2[,1], main = 'Direct'); abline(v=0.611, col = 'green')
hist(res2[,2], main = 'Indirect'); abline(v=1.26, col='green')
hist(res2[,3], main = 'Total'); abline(v=1.87, col='green')
par(og)

# SEE THE MISSING TERM 
plot(res[,3], res2[,3]); abline(0,1)
