#' Draw samples from the posterior predictive distribution
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
#' @param car_parts Data for CAR model specification; only required for \code{\link[geostan]{stan_car}} with `family = auto_gaussian()`. 
#' 
#' @param seed A single integer value to be used in a call to \code{\link[base]{set.seed}} before taking samples from the posterior distribution. 
#' 
#' @return A matrix of size S x N containing samples from the posterior predictive distribution, where S is the number of samples drawn and N is the number of observations. If `summary = TRUE`, a `data.frame` with N rows and 3 columns is returned (with column names `mu`, `lwr`, and `upr`).
#'
#' 
#' @examples
#'  fit <- stan_glm(sents ~ offset(log(expected_sents)),
#'                   re = ~ name,
#'                   data = sentencing,
#'                   family = poisson(),
#'                   chains = 2, iter = 600) # for speed only
#' 
#'  yrep <- posterior_predict(fit, S = 65)
#'  plot(density(yrep[1,]))
#'  for (i in 2:nrow(yrep)) lines(density(yrep[i,]), col = 'gray30')
#'  lines(density(sentencing$sents), col = 'darkred', lwd = 2)
#' @export
posterior_predict <- function(object, 
                S,
                summary = FALSE,
                width = 0.95,
                car_parts,
                seed
                ) {
    stopifnot(inherits(object, "geostan_fit"))
    family <- object$family$family    
    if (!missing(seed)) set.seed(seed)
    mu <- as.matrix(object, pars = "fitted")
    M <- nrow(mu)                                      
    if (missing(S)) S <- M
    if (S > M) {
        warning (paste0("Cannot draw more samples than were taken from the posterior. Using S = ", M))
        S <- M
    }
    idx <- sample(M, S)
    mu <- as.matrix(object, pars = "fitted")[idx,] 
    if (family == "auto_gaussian") {
        stopifnot(!missing(car_parts))
        rho <- as.matrix(object, "car_rho")[idx, ]
        tau <- as.matrix(object, "car_scale")[idx, ]
        preds <- .pp_auto_gaussian(mu, rho, tau, car_parts)
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

#' @importFrom Matrix Diagonal Matrix
.pp_auto_gaussian <- function(mu, rho, tau, car_parts) {
    I <- Matrix::Diagonal(car_parts$n)
    M <- Matrix::Diagonal(x = 1 / car_parts$Delta_inv)
    C <- Matrix::Matrix(car_parts$C)
    t(sapply(1:nrow(mu), function(s) {
        Sigma <- solve(I - rho[s] * C) %*% M * tau[s]^2
        MASS::mvrnorm(n = 1, mu = mu[s,], Sigma = Sigma)
    }))        
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


