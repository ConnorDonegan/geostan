#' Sample from the posterior predictive distribution
#'
#' @description Draw samples from the posterior predictive distribution of a fitted \code{geostan} model. 
#' 
#' @param object A \code{geostan_fit} object.
#' 
#' @param S Optional; number of samples to take from the posterior distribution. The default, and maximum, is the total number of samples stored in the model.
#' 
#' @param summary Should the predictive distribution be summarized by its means and central quantile intervals? If \code{summary = FALSE}, an `S` x `N` matrix of samples will be returned. If \code{summary = TRUE}, then a `data.frame` with the means and `100*width` credible intervals is returned.
#' 
#' @param width Only used if \code{summary = TRUE}, to set the quantiles for the credible intervals. Defaults to `width = 0.95`.
#'
#' @param approx For SAR models only; `approx = TRUE` uses an approximation method for the inverse of matrix `(I - rho * W)`.
#'
#' @param K For SAR models only; number of matrix powers to for the matrix inverse approximation (used when `approx = TRUE`). High values of rho (especially > 0.9) require larger K for accurate approximation.
#'
#' @param preserve_order If `TRUE`, the order of posterior draws will remain fixed; the default is to permute the MCMC samples so that (with small sample size `S`) each successive call to `posterior_predict` will return a different sample from the posterior probability distribution. 
#' 
#' @param seed A single integer value to be used in a call to \code{\link[base]{set.seed}} before taking samples from the posterior distribution. 
#' 
#' @return A matrix of size S x N containing samples from the posterior predictive distribution, where S is the number of samples drawn and N is the number of observations. If `summary = TRUE`, a `data.frame` with N rows and 3 columns is returned (with column names `mu`, `lwr`, and `upr`).
#'
#' @details
#'
#' This method returns samples from the posterior predictive distribution of the model (at the observed values of covariates, etc.). The predictions incorporate uncertainty of all parameter values (used to calculate the expected value of the model, for example) plus the error term (the model's description of the amount of variability of observations around the expected value). If the model includes measurement error in the covariates, this source of uncertainty (about \eqn{X}) is passed into the posterior predictive distribution as well.
#' 
#' For SAR models (and all other models), the observed outcomes are *not* used to formulate the posterior predictive distribution. The posterior predictive distribution for the SLM (see \link[geostan]{stan_sar}) is given by
#' \deqn{(I - \rho W)^{-1} (\mu + \epsilon).}
#' The SDLM is the same but includes spatially-lagged covariates in \eqn{mu}. The `approx = FALSE` method for SAR models requires a call to `Matrix::solve(I - rho * W)` for each MCMC sample; the `approx = TRUE` method uses an approximation based on matrix powers (LeSage and Pace 2009). The approximation will deteriorate if \eqn{\rho^K} is not near zero, so use with care.
#'
#' @source
#' LeSage, James, & Robert kelley Pace (2009). *Introduction to Spatial Econometrics*. Chapman and Hall/CRC.
#'
#' Gelman, A., J. B.Carlin, H. S. Stern, D. B. Dunson, A. Vehtari, & D. B. Rubin, D. B. (2014). *Bayesian data analysis* (3rd ed.). CRC Press.
#'
#' McElreath, Richard (2016). *Statistical Rethinking: A Bayesian Course with Examples in R and Stan*. CRC Press, Ch. 3.
#' 
#' @examples
#' E <- sentencing$expected_sents
#' sentencing$log_E <- log(E)
#'  fit <- stan_glm(sents ~ offset(log_E),
#'                   re = ~ name,
#'                   data = sentencing,
#'                   family = poisson(),
#'                   chains = 2, iter = 600) # for speed only
#'
#' 
#'  yrep <- posterior_predict(fit, S = 65)
#'  plot(density(yrep[1,] / E ))
#'  for (i in 2:nrow(yrep)) lines(density(yrep[i,] / E), col = 'gray30')
#'  lines(density(sentencing$sents / E), col = 'darkred', lwd = 2)
#'
#' sars <- prep_sar_data2(row = 9, col = 9)
#' W <- sars$W
#' y <- sim_sar(rho = .9, w = W)
#' fit <- stan_sar(y ~ 1, data = data.frame(y=y), sar = sars,
#'                 iter = 650, quiet = TRUE)
#' yrep <- posterior_predict(fit, S = 15)
#' 
#' @export
posterior_predict <- function(object, 
                              S,
                              summary = FALSE,
                              width = 0.95,
                              approx = TRUE,
                              K = 20,
                              preserve_order = FALSE,
                              seed
                              ) {
    stopifnot(inherits(object, "geostan_fit"))
    family <- object$family$family    
    if (!missing(seed)) set.seed(seed)

    mu <- as.matrix(object, pars = "fitted")    
    M <- nrow(mu)
    
    if (missing(S)) {
        S <- M
    }
    
    if (S > M) {
        warning (paste0("Cannot draw more samples than were taken from the posterior. Using S = ", M))
            S <- M
    }

    if (preserve_order) {
        idx <- seq(S)
    } else {
        idx <- sample(M, size = S)
    }

    mu <- mu[idx,]
    
    if (family == "auto_gaussian") {

        if (object$spatial$method == "CAR") {            
            rho <- as.matrix(object, "car_rho")[idx, ]
            tau <- as.matrix(object, "car_scale")[idx, ]
            preds <- .pp_car_normal(mu, rho, tau, object$car_parts)
        }

        if (object$spatial$method == "SAR") {            
            rho <- as.matrix(object, "sar_rho")[idx, ]
            tau <- as.matrix(object, "sar_scale")[idx, ]
            preds <- .pp_sar_normal(mu, rho, tau, object$sar_parts$W, approx = approx, K = K, object$sar_type)
        }
        
    }
    
    if (family == "gaussian") {
        sigma <- as.matrix(object, pars = "sigma")[idx,]
        preds <- .pp_gaussian(mu, sigma)
    }
    
    if (family == "student_t") {
        sigma <- as.matrix(object, pars = "sigma")[idx,]
        nu <- as.matrix(object, pars = "nu")[idx,]
        preds <- .pp_student(nu, mu, sigma)
    }
    
    if (family == "poisson") preds <- .pp_poisson(mu)
    
    if (family == "binomial") {
        trials <- object$data[,1] + object$data[,2]
        preds <- .pp_binomial(mu, trials)
    }
    
    if (summary) {
        df <- .pp_summary(preds, width = width)
        return(df)
    }
    
    return(preds)    
}


#' @importFrom Matrix Diagonal Matrix solve
.pp_car_normal <- function(mu, rho, tau, car_parts) {
    I <- Matrix::Diagonal(car_parts$n)
    M <- Matrix::Diagonal(x = 1 / car_parts$Delta_inv)
    C <- Matrix::Matrix(car_parts$C)
    t(sapply(1:nrow(mu), function(s) {
        Sigma <- Matrix::solve(I - rho[ s ] * C) %*% M * tau[ s ]^2
        MASS::mvrnorm(n = 1, mu = mu[ s, ], Sigma = Sigma)
    }))        
}


.pp_sar_normal <- function(mu, rho, sigma, W, approx, K, type) {

    if (!approx) {
        Draws <- sapply(1:nrow(mu), function(s) {
            sim_sar(mu =     mu[ s, ],
                    rho =   rho[ s  ],
                    sigma = sigma[ s  ],
                    w = W,
                    type = type,
                    approx = FALSE)
        }) |>
            t()
        
    return (Draws)
    }
    
    N <- nrow(W)
    S <- nrow(mu)
    Q <- K + 1
    P = 0:K
    I <- Matrix::diag(N)
    M <- list()
    M[[1]] <- I
    M[[2]] <- W
    W_k <- W    
    for (q in 3:Q) {
        W_k <- W %*% W_k
        M[[q]] <- W_k
    }        

    .rsem <- function(s) {
        eps <- rnorm(N, mean = 0, sd = sigma[ S ])
        rho_powers <- rho[ s ]^P
        Mpowers <- lapply(seq(Q), function(j) M[[j]] * rho_powers[j])
        Multiplier <- Reduce(`+`, Mpowers)
        res <- mu[ s, ] + (Multiplier %*% eps)[,1]
        return (res)
    }

    .rslm <- function(s) {
        eps <- rnorm(N, mean = 0, sd = sigma[ S ])
        rho_powers <- rho[ s ]^P
        Mpowers <- lapply(seq(Q), function(j) M[[j]] * rho_powers[j])
        Multiplier <- Reduce(`+`, Mpowers)
        z <- mu[ s, ] + eps
        res <- (Multiplier %*% z)[,1]
        return (res)
    }

    SLM <- grepl("SLM|SDLM", type)

    if (SLM) {
        
        Draws <- sapply(1:S, .rslm) |>
            t()
        
    } else {
        
        Draws <- sapply(1:S, .rsem) |>
            t()
        
    }    
    
    return (Draws)
}


.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}

.pp_student <- function(nu, mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
   mu[s,] + rt(ncol(mu), df = nu[s]) * sigma[s]
  }))
}

.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s, ])
  }))
}

.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s, ]) 
  }))
}


.pp_summary <- function(samples, width) {
    p.lwr <- (1 - width) / 2
    probs <- c(p.lwr, 1 - p.lwr)
    mu <- apply(samples, 2, mean)
    int <- apply(samples, 2, quantile, probs = probs)
    df <- data.frame(
        mu = mu,
        lwr = int[1,],
        upr = int[2,]
    )
    attributes(df)$interval <- list(width = width, probs = probs)
    return(df)
}


