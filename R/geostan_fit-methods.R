#' geostan_fit methods
#'
#' @description Methods for fitted geostan models: extract residuals, fitted values, posterior predictive distribution or spatial component from a spatial regression model; extract samples from the posterior distribution; print regression results; plot posterior distributions.
#'
#' @param object A fitted model object of class \code{geostan_fit}.
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param summary Logical; should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples from the posterior distribution at each observation is returned.
#' @param newdata A data frame in which to look for variables with which to predict, presumably for the purpose of viewing marginal effects. Note that if the model formula includes an offset term, `newdata` must contain a column with the appropriate name for the offset, even though the values will be ignored (you may set all values to 1); you must use the `alpha` argument to include any additional terms. Note also that any spatially-lagged covariate terms will be ignored if they were provided using the `slx` argument. If covariates in the model were centered using the `centerx` argument, the `predict.geostan_fit` method will automatically center the predictors in `newdata` internally using the values stored in `fit$x_center`. If `newdata` is missing, user arguments will be passed to the `fitted.geostan_fit` method to return the fitted values of the model.
#' 
#' @param alpha A single numeric value or a numeric vector with length equal to `nrow(newdata)`; `alpha` serves as the intercept in the linear predictor. The default is to use the posterior mean of the intercept. Even if \code{type = "response"}, this needs to be provided on the scale of the linear predictor. See `Details` for additional information.
#'
#' @param center May be a vector of numeric values or a logical scalar to pass to \code{\link[base]{scale}}. Defaults to using `object$x_center`. If the model was fit using `centerx = TRUE`, then covariates were centered and their mean values are stored in `object$x_center` and the `predict` method will use them to automatically center `newdata`; if the model was fit with `centerx = FALSE`, then `object$x_center = FALSE` and `newdata` will not be centered. 
#' 
#' @param type By default, results from `predict` are on the scale of the linear predictor (`type = "link")`). The alternative (`type = "response"`) is on the scale of the response variable. For example, the default return values for a Poisson model are log-rates, and using `type = "response"` will return the rates (by exponentiating the log-rates).
#' 
#' 
#' @param pars parameters to include; a character string (or vector) of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include histograms ("hist"), MCMC traceplots ("trace"), and density plots ("dens"). Diagnostic plots are also available such as Rhat statistics ("rhat"), effective sample size ("ess"), and MCMC autocorrelation ("ac").
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' @param fill fill color for histograms and density plots.
#' @param rates For Poisson and Binomial models, should the fitted values be returned as rates, as opposed to raw counts? Defaults to `TRUE`.
#' @param detrend For CAR models with Gaussian likelihood only (auto-gaussian); if `detrend = TRUE`, the implicit spatial trend will be removed from the residuals. The implicit spatial trend is `Trend = rho * C %*% (Y - Mu)` (see \code{\link[geostan]{stan_car}}). I.e., `resid = Y - (Mu + Trend)`.
#' @param ... additional arguments.
#'
#' @details
#' 
#' ### predict.geostan_fit
#'
#' The purpose of the predict method is to explore marginal effects of (combinations of) covariates. The method sets the intercept equal to its posterior mean (i.e., `alpha = mean(as.matrix(object, pars = "intercept"))`); the only source of uncertainty in the results is the posterior distribution of the coefficients, which can be obtained using `Beta = as.matrix(object, pars = "beta")`. The results returned by `predict.geostan_fit` are obtain by (a summary of):
#'```
#'   for (m in 1:M) preds[m,] = alpha + X * Beta[m,] 
#'```
#' where `M` is the number of MCMC samples in the model (`M = nrow(Beta)`) and `preds` is a matrix of predicted values.
#' 
#' Be aware that in non-linear models (including Poisson and Binomial models) marginal effects of each covariate are sensitive to the level of other covariates in the model. If the model includes any spatially-lagged covariates (introduced using the `slx` argument) or a spatial autocorrelation term, these terms will essentially be fixed at zero for the purposes of calculating marginal effects. To explore the impact of these (missing) terms, you can add their values to the linear predictor using the `alpha` argument. 
#' 
#' @return
#'
#' Methods \code{residuals}, \code{fitted}, \code{predict}, and \code{spatial} return a matrix containing all samples for each observation if \code{summary = FALSE}, else if \code{summary = TRUE} a \code{data.frame} containing a summary of the posterior distribution at each observation (of, respectively, residuals, fitted values, predicted values, or the spatial trend). The \code{predict} method will return a data frame with a summary of results together with use-provided `newdata`.
#'
#' The \code{predict} method is designed for reviewing marginal effects of covariates. Thus, results do not include spatial trends or offset terms. To obtain the fitted values of the model (as opposed to predictions from new data), use the \code{fitted} method. For the posterior predictive distribution, see \code{\link[geostan]{posterior_predict}}. 
#'
#' \code{plot} returns a \code{ggplot} object that can be customized using the \code{ggplot2} package.
#'
#' \code{as.matrix}, \code{as.data.frame}, \code{as.array} return samples from the joint posterior distribution of parameters in the format corresponding to their names. The \code{pars} argument is used to return samples from only a subset of parameters.
#'
#' @seealso \code{\link[geostan]{posterior_predict}}, \code{\link[geostan]{stan_glm}}, \code{\link[geostan]{stan_esf}}, \code{\link[geostan]{stan_icar}}, \code{\link[geostan]{stan_car}}
#' 
#' @examples 
#' \donttest{
#' library(ggplot2)
#' library(sf)
#' data(sentencing)
#'
#' # spatial weights matrix with binary coding scheme
#' C <- shape2mat(sentencing, style = "B")
#'
#' # log-expected number of sentences
#' ## expected counts are based on county racial composition and mean sentencing rates
#' log_e <- log(sentencing$expected_sents)
#'
#' # fit spatial Poisson model with unstructured 'random effects'
#' fit <- stan_glm(sents ~ offset(log_e),
#'                    re = ~ name,
#'                    family = poisson(),
#'                    data = sentencing,
#'                    C = C,
#'                    chains = 2, iter = 500) # for speed only
#'
#' # print and plot results
#' print(fit)
#' plot(fit)
#'
#' # residuals
#' r = resid(fit)
#'
#' # fitted values
#' f = fitted(fit)
#'
#' # spatial diagnostics
#' sp_diag(fit, sentencing)
#' 
#' # county `random effects' 
#' sp = spatial(fit)
#'
#' # posterior predictive distribution
#' yrep <- posterior_predict(fit, S = 100)
#' bayesplot::ppc_dens_overlay(sentencing$sents, yrep)
#'
#' # extract matrix of samples from posterior distribution of parameters
#' ## alpha_re are the unstructured area random effects
#' S.matrix <- as.matrix(fit, pars = "alpha_re")
#'
#' # array of samples
#' S.array <- as.array(fit, pars = c("intercept", "alpha_re", "alpha_tau"))
#' S.monitor <- rstan::monitor(S.array, print = FALSE, warmup = 0)
#' head(S.monitor)
#'
#' ## marginal effects
#' data(georgia)
#' C <- shape2mat(georgia, style = "B")
#' cp <- prep_car_data(C)
#' georgia$income <- georgia$income/1e3
#' 
#' fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)) + log(income),
#'                 slx = ~ log(income),
#'                 centerx = TRUE,
#'                 car_parts = cp,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 2, iter = 500) # for speed only
#'
#' newdata <- data.frame(
#'     income = seq(min(georgia$income), max(georgia$income), by = 1),
#'     pop.at.risk.male = 1
#' )
#'
#' p <- predict(fit, newdata, type = "response")
#' plot(newdata$income, p$mean * 1e3,
#'      main = "Deaths per 1,000",
#'      ylab = NA,
#'      xlab = "Median county income ($1,000s)")
#'
#' }
#' 
#' @export
#' @md
#' @method print geostan_fit
#' @rdname geostan_fit
print.geostan_fit <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3, pars = NULL, ...) {
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
  if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])   
  if(inherits(x$re$formula, "formula")) {
    cat("Partial pooling (varying intercept): ")
    print(x$re$formula)
    pars <- c(pars, "alpha_tau")
  }
  cat("Spatial method (outcome): ", as.character(x$spatial$method), "\n")
  if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
  if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
  if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
  if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")  
  cat("Likelihood function: ", x$family$family, "\n")
  cat("Link function: ", x$family$link, "\n")
  cat("Residual Moran Coefficient: ", x$diagnostic["Residual_MC"], "\n")
  if (!is.na(x$diagnostic["WAIC"])) cat("WAIC: ", x$diagnostic["WAIC"], "\n")
  cat("Observations: ", nrow(x$data), "\n")  
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
  if (x$spatial$method == "ESF") {
    cat("\nHorseshoe global shrinkage prior: ", round(x$priors$beta_ev$global_scale, 2), "\n")
  }
  cat("\n")
  print(x$stanfit, pars = pars, digits = digits, probs = probs, ...)
}


#' @export
#' @import graphics
#' @importFrom signs signs
#' @rdname geostan_fit
#' @method plot geostan_fit
plot.geostan_fit <- function(x, pars, plotfun = "hist", fill = "steelblue4", ...) {
  if(missing(pars)) {
      pars <- "intercept"
      if (inherits(x$slx, "formula")) pars <- c(pars, "gamma")      
      x.pars <- c("beta", "nu", "sigma", "alpha_tau")
      if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])
      if (x$spatial$method == "CAR") pars <- c(pars, "car_rho", "car_scale")
      if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
      if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
      if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")
  }
  rstan::plot(x$stanfit, pars = pars, fill = fill, plotfun = plotfun, ...) + scale_x_continuous(labels = signs::signs) + scale_y_continuous(labels = signs::signs)
}


#' @export
#' @rdname geostan_fit
#' @method as.matrix geostan_fit
as.matrix.geostan_fit <- function(x, ...) {
  as.matrix(x$stanfit, ...)
}

#' @export
#' @rdname geostan_fit
#' @method as.data.frame geostan_fit
as.data.frame.geostan_fit <- function(x, ...) {
    as.data.frame(x$stanfit, ...)
}

#' @export
#' @rdname geostan_fit
#' @method as.array geostan_fit
as.array.geostan_fit <- function(x, ...){
    as.array(x$stanfit, ...)
}


#' @noRd
.resid <- function(yhat, y) y - yhat


#' @method residuals geostan_fit
#' @rdname geostan_fit
#' @export
residuals.geostan_fit <- function(object, summary = TRUE, rates = TRUE, detrend = TRUE, ...) {
    y <- object$data[,1]
    fits <- fitted(object, summary = FALSE, rates = rates)
    if (rates && object$family$family == "binomial") { y <- y / (y + object$data[,2]) }
    if (rates && object$family$family == "poisson" && "offset" %in% c(colnames(object$data))) {
        log.at.risk <- object$data[, "offset"]
        at.risk <- exp( log.at.risk )
        y <- y / at.risk     
    }
    R = sweep(fits, MARGIN = 2, STATS = as.array(y), FUN = .resid)
    if (object$family$family == "auto_gaussian" && detrend) R <- R - spatial(object, summary = FALSE)
    colnames(R) <- gsub("fitted", "residual", colnames(R))
    if (summary) return( post_summary(R) )
    return (R)  
}

#' @export
#' 
#' @import stats
#' @method fitted geostan_fit
#' @rdname geostan_fit
fitted.geostan_fit <- function(object, summary = TRUE, rates = TRUE, ...) {
    fits <- as.matrix(object$stanfit, pars = "fitted")
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
#' @rdname geostan_fit
spatial <- function(object, summary = TRUE, ...) {
    UseMethod("spatial", object)
}

#' @export
#' @method spatial geostan_fit
#' @rdname geostan_fit
spatial.geostan_fit <- function(object, summary = TRUE, ...) {
  if (is.na(object$spatial$par)) stop("This model does not have a spatial trend component to extract.")
  par <- as.character(object$spatial$par)
  if (!object$spatial$method == "CAR") {
      spatial.samples <- as.matrix(object, pars = par, ...)
  } else {
      if (object$family$family == "auto_gaussian") {
          C <- object$C              
          R <- resid(object, summary = FALSE, detrend = FALSE)
          rho <- as.matrix(object, pars = "car_rho")
          spatial.samples <- t(sapply(1:nrow(rho), function(i) {
              as.numeric( rho[i] * C %*% R[i,] )
          }))         
      } else {
          log_lambda_mu <- as.matrix(object, pars = "log_lambda_mu")
          log_lambda <- log( fitted(object, summary = FALSE, rates = TRUE) )
          spatial.samples <- log_lambda - log_lambda_mu
      }
  } 
  if (summary) {
    return ( post_summary(spatial.samples) )
  } else {
      return( spatial.samples )
  }
}

#' @export
#' @method predict geostan_fit
#' @rdname geostan_fit
predict.geostan_fit <- function(object, newdata, alpha = mean(as.matrix(object, pars = "intercept")), center = object$x_center, summary = TRUE, type = c("link", "response"), ...) {
    type <- match.arg(type)
    if (missing(newdata)) return (fitted(object, summary = summary, ...))
    f <- object$formula[-2]
    X <- as.matrix(model.matrix(f, newdata)[,-1])    
    X <- scale(X, center = center, scale = FALSE)
    B <- as.matrix(object, pars = "beta") 
    M <- nrow(B)
    N <- nrow(X)
    P <- matrix(NA, nrow = M, ncol = N)
    for (m in 1:M) P[m,] <- alpha + X %*% B[m,] 
    if (type == "response") {
        if (object$family$link == "log") P <- exp(P)
        if (object$family$link == "logit") P <- inv_logit(P)
    }
    if (summary) {
        P <- post_summary(P)
        P <- cbind(newdata, alpha = alpha, P)
    }
    return (P)
}

