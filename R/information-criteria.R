
#' Model comparison
#'
#' @description Deviance Information Criteria (DIC) and Widely Application Information Criteria (WAIC) for model comparison.
#' 
#' @param object A fitted \code{geostan} model
#' 
#' @param pointwise Logical (defaults to `FALSE`), should a vector of values for each observation be returned?
#' 
#' @param digits Round results to this many digits.
#' 
#' @return
#'
#' WAIC returns a vector of length 3 with the \code{WAIC} value, a penalty term which measures the effective number of parameters estimated by the model \code{Eff_pars}, and log predictive density \code{Lpd}. If \code{pointwise = TRUE}, results are returned in a \code{data.frame}.
#'
#' DIC returns a vector of length 2: the DIC value and the penalty term (which is part of the DIC calculation).
#'
#' @details
#' 
#' WAIC (widely applicable information criteria) and DIC (deviance information criteria) are used for model comparison. They are based on theories of out-of-sample predictive accuracy. The DIC is implemented with penalty term defined as 1/2 times the posterior variance of the deviance (Spiegelhatler et al. 2014).
#'
#' The limitations of these methods include that DIC is less robust than WAIC and that WAIC is not strictly valid for autocorrelated data (viz. geostan's spatial models).
#'
#' For both DIC and WAIC, lower values indicate better models. 
#' 
#' @examples
#' data(georgia)
#' 
#' fit <- stan_glm(log(rate.male) ~ 1, data = georgia,
#'                 iter=600, chains = 2, quiet = TRUE)
#' fit2 <- stan_glm(log(rate.male) ~ log(income), data = georgia,
#'                 centerx = TRUE, iter=600, chains = 2, quiet = TRUE)
#' 
#' dic(fit)
#' dic(fit2)
#' 
#' waic(fit)
#' waic(fit2)
#' @source
#'
#' D. Spiegelhatler, N. G. Best, B. P. Carlin and G. Linde (2014) The Deviance Information Criterion: 12 Years on. J. Royal Statistical Society Series B: Stat Methodology. 76(3): 485-493.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information criterion in singular learning theory. Journal of Machine Learning Research 11, 3571-3594.
#'
#' @md 
#' @export
#' @rdname waic
waic <- function(object, pointwise = FALSE, digits = 2) {
  ll <- log_lik(object, array = FALSE)
  nsamples <- nrow(ll)
  lpd <- apply(ll, 2, log_sum_exp) - log(nsamples)
  p_waic <- apply(ll, 2, var)
  waic <- -2 * (lpd - p_waic)
  if(pointwise) return(data.frame(waic = waic, eff_pars = p_waic, lpd = lpd))
  res <- c(WAIC = sum(waic), Eff_pars = sum(p_waic), Lpd = sum(lpd))
  return(round(res, digits))
}


#' @export
#' @rdname waic
dic <- function(object, digits = 1) {
    
    # get log-likelihood
    ll <- log_lik(object, array = FALSE)

    # calculate model deviance
    dev <- -2 * apply(ll,  1, FUN = sum) 
    dev_bar <- mean(dev)
    
    # calculate penalty term
    penalty <- 0.5 * var(dev)
    
    # DIC
    DIC <- dev_bar + penalty

    # round digits
    x = round(c(DIC = DIC, penalty = penalty), digits)

    # return DIC and penalty term
    return (x)
}

#' Extract log-likelihood
#'
#' @param object A `geostan_fit` model
#' 
#' @param array Return results as an array, one matrix per MCMC chain?
#' 
#' @param ... Other arguments (not used)
#'
#' @return A matrix (or array) of MCMC samples for the log-likelihood: the casewise probability of the data conditional on estimated parameter values.
#'
#' @seealso \code{\link[geostan]{waic}} \code{\link[geostan]{dic}}
#' 
#' @export
#' @rdname log_lik
log_lik <- function(object, array = FALSE, ...) {
        UseMethod("log_lik", object)
    }


#' @export
#' @method log_lik geostan_fit
#' @importFrom stats dt dnorm dpois
#' @rdname log_lik
log_lik.geostan_fit <- function(object, array = FALSE, ...) {
    
    M <- object$stanfit@sim$iter - object$stanfit@sim$warmup
    K <- object$stanfit@sim$chains
    N <- object$N
    
    mu <- as.array(object, pars = "fitted")
    idx <- object$missing$y_obs_idx
    fam <- object$family$family
    
    if (fam == "student_t") {
        log_lik <- array(dim = c(M, K, N))
        y <- object$data[,'y']        
        sigma <- as.array(object, pars = "sigma")
        nu <- as.array(object, pars = "nu")        
        for (m in 1:M) {
            for (k in 1:K) {                    
                log_lik[m,k,] <- stats::dt(
                    x = (y[idx] - mu[m,k,idx])/sigma[m,k,1],
                    df = nu[m,k,1],
                    log = TRUE) - log(sigma[m,k,1])                
            }
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)            
    }
    
    if (fam == "gaussian") {
        log_lik <- array(dim = c(M, K, N))        
        y <- object$data[,'y']        
        sigma <- as.array(object, pars = "sigma")
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,] <- dnorm(
                    x = y[idx],
                    mean = mu[m,k,idx],
                    sd = sigma[m,k,1],
                    log = TRUE)
            }
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)           
    }
    
    if (fam == "poisson") {
        cp <- object$missing$censor_point
        if (cp == FALSE) {            
            log_lik <- array(dim = c(M, K, N))
            y <- object$data[,'y']        
            for (m in 1:M) {
                for (k in 1:K) {                    
                    log_lik[m,k,] <- dpois(
                        x = y[idx],
                        lambda = mu[m,k,idx],
                        log = TRUE)
                }
            }                               
        } else {
            N <- object$missing$n_obs + object$missing$n_mis
            log_lik <- array(dim = c(M, K, N))
            y <- object$data[,'y']
            mis_idx <- object$missing$y_mis_idx
            for (m in 1:M) {
                for (k in 1:K) {                    
                    log_lik[m,k,idx] <- dpois(
                        x = y[idx],
                        lambda = mu[m,k,idx],
                        log = TRUE)
                    log_lik[m,k,mis_idx] <- ppois(q = cp,
                                                  lambda = mu[m,k,mis_idx],
                                                  log.p = TRUE)
                }
            }            
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)            
    }
    
    if (fam == "binomial") {
        log_lik <- array(dim = c(M, K, N))
        #out = model.response(model.frame(object$formula, object$data, na.action = NULL))
        dat <- object$data
        y <- dat[,1]
        size <- y + dat[,2]        
        for (m in 1:M) {
            for (k in 1:K) {                    
                log_lik[m,k,] <- dbinom(
                    x = y[idx],
                    size = size[idx],
                    prob = mu[m,k,idx],
                    log = TRUE)
            }
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)            
    }

    if (fam == "auto_gaussian" & object$spatial$method == "CAR") {
        log_lik <- array(dim = c(M, K, 1))
        y <- object$data[,'y']
        sigma <- as.array(object, pars = "car_scale")
        rho <- as.array(object, pars = "car_rho")
        parts <- object$car_parts
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,1] <- car_normal_lpdf(
                    y = y,
                    mu = mu[m,k,],
                    sigma = sigma[m,k,1],
                    rho = rho[m,k,1],
                    C = parts$C,
                    D_inv = parts$Delta_inv,
                    log_det_D_inv = parts$log_det_Delta_inv,
                    lambda = parts$lambda,
                    n = parts$n
                )
            }
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = 1)            
    }

    if (fam == "auto_gaussian" & object$spatial$method == "SAR") {
        log_lik <- array(dim = c(M, K, 1))
        y <- object$data[,'y']
        sigma <- as.array(object, pars = "sar_scale")
        rho <- as.array(object, pars = "sar_rho")
        parts <- object$sar_parts
        for (m in 1:M) {
            for (k in 1:K) {
                log_lik[m,k,1] <- sar_normal_lpdf(
                    y = y,
                    mu = mu[m,k,],
                    sigma = sigma[m,k,1],
                    rho = rho[m,k,1],
                    W = parts$W,
                    lambda = parts$eigenvalues_w,
                    n = parts$n
                )
            }
        }
        if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = 1)            
    }
    
    return (log_lik)
}
