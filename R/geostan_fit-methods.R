#' print or plot a fitted geostan model
#'
#' @description Print a summary of model results to the R console, or plot posterior distributions of model parameters.
#'
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' 
#' @param pars parameters to include; a character string (or vector) of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include histograms ("hist"), MCMC traceplots ("trace"), and density plots ("dens"). Diagnostic plots are also available such as Rhat statistics ("rhat"), effective sample size ("ess"), and MCMC autocorrelation ("ac").
#' 
#' @param fill fill color for histograms and density plots.
#'
#' @param ... additional arguments to `rstan::plot` or `rstan::print.stanfit`.
#'
#' @return
#'
#' The print methods writes text to the console to summarize model results. The plot method resturns a `ggplot` (from `rstan::plot` for stanfit objects).
#' 
#' @examples
#' data(georgia)
#' georgia$income <- georgia$income/1e3
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + log(income),
#'                 centerx = TRUE,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 2, iter = 600) # for speed only
#'
#' 
#' # print and plot results
#' print(fit)
#' plot(fit)
#' @export
#' @md
#' @method print geostan_fit
#' @rdname print_geostan_fit
print.geostan_fit <- function(x,
                              probs = c(0.025, 0.2, 0.5, 0.8, 0.975),
                              digits = 3,
                              pars = NULL, ...) {
  pars <- c(pars, "intercept")
  cat("Spatial Model Results \n")
  cat("Formula: ")
  print(x$formula)
  if (inherits(x$slx, "formula")) {
   cat("SLX: ")
   print(x$slx)
   pars <- c(pars, "gamma")
  }
  x.pars <- c("beta", "nu", "sigma")
  if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="$|"), names(x$priors))])   
  if(inherits(x$re$formula, "formula")) {
    cat("Partial pooling (varying intercept): ")
    print(x$re$formula)
    pars <- c(pars, "alpha_tau")
  }  
  cat("Spatial method (outcome): ", as.character(x$spatial$method), "\n")
  if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
  if (x$spatial$method == "SAR") pars <- c(pars, "sar_rho", "sar_scale")  
  if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
  if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
  if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")  
  cat("Likelihood function: ", x$family$family, "\n")
  cat("Link function: ", x$family$link, "\n")
  if (!is.null(x$diagnostic$Residual_MC)) cat("Residual Moran Coefficient: ", x$diagnostic$Residual_MC, "\n")
  if (!is.null(x$diagnostic$WAIC)) cat("WAIC: ", x$diagnostic$WAIC, "\n")
  cat("Observations: ", x$N, "\n")  
  cat("Data models (ME): ")
  if (x$ME$has_me) {
      cat(paste(names(x$ME$se), sep = ", ")) 
      if (x$ME$spatial_me) {
          cat("\n Data model (ME prior): CAR (auto Gaussian)")
      } else {
          cat("\nData model (ME prior): Student's t")
      }
  } else {
      cat("none")
  }
  if (x$spatial$method == "HS") {
    cat("\nHorseshoe global shrinkage prior: ", round(x$priors$beta_ev$global_scale, 2), "\n")
  }
  cat("\n")
  print(x$stanfit, pars = pars, digits = digits, probs = probs, ...)
}


#' @export
#' @import graphics
#' @importFrom signs signs
#' @rdname print_geostan_fit
#' @method plot geostan_fit
plot.geostan_fit <- function(x,
                             pars,
                             plotfun = "hist",
                             fill = "steelblue4",
                             ...) {
  if(missing(pars)) {
      pars <- "intercept"
      if (inherits(x$slx, "formula")) pars <- c(pars, "gamma")      
      x.pars <- c("beta", "nu", "sigma", "alpha_tau")
      if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])
      if (x$spatial$method == "SAR") pars <- c(pars, "sar_rho", "sar_scale")  
      if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
      if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
      if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
      if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")
  }
  rstan::plot(x$stanfit, pars = pars, fill = fill, plotfun = plotfun, ...) +
      scale_x_continuous(labels = signs::signs) +
      scale_y_continuous(labels = signs::signs)
}




#' Extract residuals, fitted values, or the spatial trend 
#'
#' @description Extract model residuals, fitted values, or spatial trend from a fitted \code{geostan_fit} model. 
#'
#' @param object A fitted model object of class \code{geostan_fit}.
#' 
#' @param summary Logical; should the values be summarized by their mean, standard deviation, and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise, a matrix containing samples from the posterior distributions is returned.
#'
#' @param rates For Poisson and Binomial models, should the fitted values be returned as rates, as opposed to raw counts? Defaults to `TRUE`; see the `Details` section for more information.
#' 
#' @param detrend For auto-normal models (CAR and SAR models with Gaussian likelihood only); if `detrend = TRUE`, the implicit spatial trend will be removed from the residuals. The implicit spatial trend is `Trend = rho * C %*% (Y - Mu)` (see \link[geostan]{stan_car} or \link[geostan]{stan_sar}). I.e., `resid = Y - (Mu + Trend)`.
#' 
#' @param trend For auto-normal models (CAR and SAR models with Gaussian likelihood only); if `trend = TRUE`, the fitted values will include the implicit spatial trend term. The implicit spatial trend is `Trend = rho * C %*% (Y - Mu)` (see \link[geostan]{stan_car} or \link[geostan]{stan_sar}). I.e., if `trend = TRUE`, `fitted = Mu + Trend`.
#' 
#' @param ... Not used
#'
#' @return
#'
#' By default, these methods return a `data.frame`. The column named `mean` is what most users will be looking for. These contain the fitted values (for the `fitted` method), the residuals (fitted values minus observed values, for the `resid` method), or the spatial trend (for the `spatial` method). The `mean` column is the posterior mean of each value, and the column `sd` contains the posterior standard deviation for each value. The posterior distributions are also summarized by select quantiles (including 2.5\% and 97.5\%). 
#'
#' If `summary = FALSE` then the method returns an S-by-N matrix of MCMC samples, where S is the number of MCMC samples and N is the number of observations in the data.
#' 
#' @details
#'
#' When \code{rates = FALSE} and the model is Poisson or Binomial, the fitted values returned by the \code{fitted} method are the expected value of the response variable. The \code{rates} argument is used to translate count outcomes to rates by dividing by the appropriate denominator. The behavior of the `rates` argument depends on the model specification. Consider a Poisson model of disease incidence, such as the following intercept-only case:
#' ```
#' fit <- stan_glm(y ~ offset(log(E)),
#'                data = data,
#'                family = poisson())
#' ```
#' If the fitted values are extracted using `rates = FALSE`, then \code{fitted(fit)} will return the expectation of \eqn{y}. If `rates = TRUE` (the default), then \code{fitted(fit)} will return the expected value of the rate \eqn{\frac{y}{E}}.
#'
#' If a binomial model is used instead of the Poisson, then using `rates = TRUE` will return the expectation of \eqn{\frac{y}{N}} where \eqn{N} is the sum of the number of 'successes' and 'failures', as in:
#' ```
#' fit <- stan_glm(cbind(successes, failures) ~ 1,
#'                data = data,
#'                family = binomial())
#' ```
#' 
#' @examples
#' \donttest{
#' data(georgia)
#' A <- shape2mat(georgia, "B")
#' 
#' fit <- stan_esf(deaths.male ~ offset(log(pop.at.risk.male)),
#'                 C = A,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 1, iter = 600) # for speed only
#'
#' 
#' # Residuals
#' r <- resid(fit)
#' moran_plot(r$mean, A)
#' head(r)
#' 
#' # Fitted values
#' f <- fitted(fit)
#'
#' # Fitted values, unstandardized
#' f <- fitted(fit, rates = FALSE)
#' head(f)
#'
#' # Spatial trend
#' esf <- spatial(fit)
#' head(esf)
#' }
#' @export
#' @md
#' @method residuals geostan_fit
#' @rdname resid_geostan_fit
#' @export
residuals.geostan_fit <- function(object,
                                  summary = TRUE,
                                  rates = TRUE,
                                  detrend = TRUE,
                                  ...) {
    y <- object$data[,1]
    fits <- fitted(object, summary = FALSE, rates = rates, trend = detrend)
    if (rates && object$family$family == "binomial") { y <- y / (y + object$data[,2]) }
    if (rates && object$family$family == "poisson" && "offset" %in% c(colnames(object$data))) {
        log.at.risk <- object$data[, "offset"]
        at.risk <- exp( log.at.risk )
        y <- y / at.risk     
    }
    R = sweep(fits, MARGIN = 2, STATS = as.array(y), FUN = .resid)
    colnames(R) <- gsub("fitted", "residual", colnames(R))
    if (summary) return( post_summary(R) )
    return (R)  
}


#' @noRd
.resid <- function(yhat, y) y - yhat

#' @export
#' 
#' @import stats
#' @method fitted geostan_fit
#' @rdname resid_geostan_fit
fitted.geostan_fit <- function(object,
                               summary = TRUE,
                               rates = TRUE,
                               trend = TRUE,
                               ...) {
    fits <- as.matrix(object$stanfit, pars = "fitted")
    if (object$family$family == "auto_gaussian" && trend == TRUE) {
        spatial_samples <- spatial(object, summary = FALSE)
        fits <- fits + spatial_samples
    }
    if (!rates && object$family$family == "binomial") {        
        trials <- object$data[,1] + object$data[,2]
        fits <- sweep(fits, MARGIN = 2, STATS = as.array(trials), FUN = "*")
    }
    if (rates && object$family$family == "poisson" && "offset" %in% c(colnames(object$data))) {
        log.at.risk <- object$data[, "offset"]
        at.risk <- exp( log.at.risk )
        fits <- sweep(fits, MARGIN = 2, STATS = as.array(at.risk), FUN = "/")
    }
    if (summary) return( post_summary(fits) )
    return (fits)          
}

#' Extract spatial component from a fitted geostan model
#' @export
#' @rdname resid_geostan_fit
spatial <- function(object,
                    summary = TRUE,
                    ...) {
    UseMethod("spatial", object)
}

#' @export
#' @method spatial geostan_fit
#' @rdname resid_geostan_fit
spatial.geostan_fit <- function(object,
                                summary = TRUE,
                                ...) {
    if (is.na(object$spatial$par)) stop("This model does not have a spatial trend component to extract.")
    if (object$spatial$method %in% c("CAR", "SAR")) {
        spatial.samples <- extract_autoGauss_trend(object)
    } else {
        par <- as.character(object$spatial$par)
        spatial.samples <- as.matrix(object, pars = par, ...)
    }
  if (summary) {
    return ( post_summary(spatial.samples) )
  } else {
      return( spatial.samples )
  }
}

extract_autoGauss_trend <- function(object) {
    if (object$spatial$method == "CAR") rho_name <- "car_rho"
    if (object$spatial$method == "SAR") rho_name <- "sar_rho"
    if (object$family$family == "auto_gaussian") {
        C <- object$C
        y <- object$data[,1]
        rho <- as.matrix(object, pars = rho_name)
        fits <- as.matrix(object$stanfit, pars = "fitted")
        R = sweep(fits, MARGIN = 2, STATS = as.array(y), FUN = .resid)
        spatial.samples <- t(sapply(1:nrow(rho), function(i) {
            as.numeric( rho[i] * C %*% R[i,] )
        }))         
      } else {
          log_lambda_mu <- as.matrix(object, pars = "log_lambda_mu")
          log_lambda <- log( fitted(object, summary = FALSE, rates = TRUE) )
          spatial.samples <- log_lambda - log_lambda_mu
      }
    return (spatial.samples)
}



#' Extract samples from a fitted model
#'
#' @description Extract samples from the joint posterior distribution of parameters.
#'
#' @param x A fitted model object of class \code{geostan_fit}.
#'
#' @param ... Further arguments passed to \code{rstan} methods for for `as.data.frame`, `as.matrix`, or `as.array`
#'
#' @return
#'
#' A matrix, data frame, or array of MCMC samples is returned.
#' 
#' @examples
#' data(georgia)
#' A <- shape2mat(georgia, "B")
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)),
#'                 C = A,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 1, iter = 600) # for speed only
#'
#' 
#' s <- as.matrix(fit)
#' dim(s)
#'
#' a <- as.matrix(fit, pars = "intercept")
#' dim(a)
#'
#' # Or extract the stanfit object
#' S <- fit$stanfit
#' print(S, pars = "intercept")
#' samples <- as.matrix(S)
#' dim(samples)
#' @md
#' @export
#' @rdname samples_geostan_fit
#' @method as.matrix geostan_fit
as.matrix.geostan_fit <- function(x, ...) {
  as.matrix(x$stanfit, ...)
}

#' @export
#' @rdname samples_geostan_fit
#' @method as.data.frame geostan_fit
as.data.frame.geostan_fit <- function(x, ...) {
    as.data.frame(x$stanfit, ...)
}

#' @export
#' @rdname samples_geostan_fit
#' @method as.array geostan_fit
as.array.geostan_fit <- function(x, ...){
    as.array(x$stanfit, ...)
}


#' Predict method for `geostan_fit` models
#'
#' @description Obtain predicted values from a fitted model by providing new covariate values.
#' 
#' @param object A fitted model object of class \code{geostan_fit}.
#' 
#' @param newdata A data frame in which to look for variables with which to predict, presumably for the purpose of viewing marginal effects. Note that if the model formula includes an offset term, `newdata` must contain the offset. Note also that any spatially-lagged covariate terms will be ignored if they were provided using the `slx` argument. If covariates in the model were centered using the `centerx` argument, the `predict.geostan_fit` method will automatically center the predictors in `newdata` using the values stored in `object$x_center`. If `newdata` is missing, the fitted values of the model will be returned.
#' 
#' @param alpha An N-by-1 matrix of MCMC samples for the intercept; this is provided by default. However, this argument might be used if there is a need to incorporate the spatial trend term, in which case it may be thought of as a spatially-varying intercept. If used, note that the intercept needs to be provided on the scale of the linear predictor. 
#'
#' @param center Optional vector of numeric values or a logical scalar to pass to \code{\link[base]{scale}}. Defaults to using `object$x_center`. If the model was fit using `centerx = TRUE`, then covariates were centered and their mean values are stored in `object$x_center` and the `predict` method will use them to automatically center `newdata`; if the model was fit with `centerx = FALSE`, then `object$x_center = FALSE` and `newdata` will not be centered.
#'
#' @param summary If `FALSE`, a matrix containing samples from the posterior distribution at each observation is returned. The default, `TRUE`, will summarize results by providing an estimate (mean) and credible interval (formed by taking quantiles of the MCMC samples).
#' 
#' @param type By default, results from `predict` are on the scale of the linear predictor (`type = "link")`). The alternative (`type = "response"`) is on the scale of the response variable. For example, the default return values for a Poisson model on the log scale, and using `type = "response"` will return the original scale of the outcome variable (by exponentiating the log values).
#'
#' @param ... Not used
#'
#' @details
#'
#' The purpose of the predict method is to explore marginal effects of (combinations of) covariates. 
#'
#' The model formula will be taken from `object$formula`, and then a model matrix will be created by passing `newdata` to the \link[stats]{model.frame} function (as in: \code{model.frame(newdata, object$formula}). Parameters are taken from `as.matrix(object, pars = c("intercept", "beta"))`.
#' 
#' Be aware that in generalized linear models (such as Poisson and Binomial models) marginal effects plots on the response scale may be sensitive to the level of other covariates in the model. If the model includes any spatially-lagged covariates (introduced using the `slx` argument) or a spatial autocorrelation term (for example, you used a spatial CAR, SAR, or ESF model), these terms will essentially be fixed at zero for the purposes of calculating marginal effects. If you want to change this, you can introduce spatial trend values by specifying a varying intercept using the `alpha` argument.
#' 
#' @return
#'
#' If `summary = FALSE`, a matrix of samples is returned. If `summary = TRUE` (the default), a data frame is returned.
#' 
#' @examples
#' data(georgia)
#' georgia$income <- georgia$income / 1e3
#' 
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + log(income),
#'                data = georgia,
#'                centerx = TRUE,
#'                family = poisson(),
#'                chains = 2, iter = 600) # for speed only
#'
#' # note: pop.at.risk.male=1 leads to log(pop.at.risk.male)=0
#' # so that the predicted values are rates
#' newdata <- data.frame(
#'              income = seq(min(georgia$income),
#'                           max(georgia$income),
#'                            length.out = 100),
#'              pop.at.risk.male = 1)
#'
#' preds <- predict(fit, newdata, type = "response")
#' head(preds)
#' plot(preds$income,
#'      preds$mean * 10e3,
#'      type = "l",
#'      ylab = "Deaths per 10,000",
#'      xlab = "Income ($1,000s)")
#'
#' # here the predictions are rates per 10,000
#' newdata$pop.at.risk.male <- 10e3
#' preds <- predict(fit, newdata, type = "response")
#' head(preds)
#' plot(preds$income,
#'     preds$mean,
#'     type = "l",
#'     ylab = "Deaths per 10,000",
#'     xlab = "Income ($1,000s)")
#' @export 
predict.geostan_fit <- function(object,
                                newdata,
                                alpha = as.matrix(object, pars = "intercept"),
                                center = object$x_center,
                                summary = TRUE,
                                type = c("link", "response"),
                                ...) {
    type <- match.arg(type)
    if (missing(newdata)) return (fitted(object, summary = summary, ...))    
    f <- object$formula[-2]
    X <- as.matrix(model.matrix(f, newdata)[,-1])    
    X <- scale(X, center = center, scale = FALSE)
    O <- model.offset( model.frame(f, newdata) )
    O <- ifelse(is.null(O), 0, O)
    Beta <- as.matrix(object, pars = "beta") 
    M <- nrow(Beta)
    N <- nrow(X)
    P <- matrix(NA, nrow = M, ncol = N)
    stopifnot(nrow(alpha) == nrow(Beta))
    for (m in 1:M) P[m,] <- O + alpha[m,] + X %*% Beta[m,] 
    if (type == "response") {
        if (object$family$link == "log") P <- exp(P)
        if (object$family$link == "logit") P <- inv_logit(P)
    }
    if (summary) {
        P <- post_summary(P)
        P <- cbind(newdata, P)
    }
    return (P)
}



#' @export
log_lik <- function(object, ...) {
        UseMethod("log_lik", object)
    }


#' @export
#' @method log_lik geostan_fit
#' @importFrom stats dt dnorm dpois
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
    }
    
    if (fam == "binomial") {
        log_lik <- array(dim = c(M, K, N))
        out = model.response(model.frame(object$formula, object$data, na.action = NULL))
        y <- out[,1]
        size <- y + out[,2]        
        for (m in 1:M) {
            for (k in 1:K) {                    
                log_lik[m,k,] <- dbinom(
                    x = y[idx],
                    size = size[idx],
                    p = mu[m,k,idx],
                    log = TRUE)
            }
        }
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
    }

    if (array == FALSE) log_lik <- matrix(log_lik, nrow = K * M, ncol = N)
    
    return (log_lik)
}


