#' geostan_fit methods
#'
#' @description Methods for fitted geostan models: extract residuals, fitted values, posterior predictive distribution or spatial component from a spatial regression model; extract samples from the posterior distribution; print regression results; plot posterior distributions.
#'
#' 
#' @param object A fitted model object of class \code{geostan_fit}.
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param summary Logical; should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples from the posterior distribution at each observation is returned.
#' @param pars parameters to include; a character string or vector of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include histograms ("hist"), MCMC traceplots ("trace"), and density plots ("dens"). Diagnostic plots are also available such as Rhat statistics ("rhat"), effective sample size ("ess"), and MCMC autocorrelation ("ac").
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' @param fill fill color for histograms and density plots.
#' @param rates For Poisson and Binomial models, should the fitted values be returned as rates, as opposed to raw counts? Defaults to `TRUE`.
#' @param detrend For CAR models with Gaussian likelihood only (auto-gaussian); if `detrend = TRUE`, the implicit spatial trend will be removed from the residuals. The implicit spatial trend is `Trend = rho * C %*% (Y - Mu)` (see \code{\link[geostan]{stan_car}}). I.e., `resid = Y - (Mu + Trend)`.
#' @param ... additional arguments.
#' @return
#'
#' Methods \code{residuals}, \code{fitted}, \code{spatial} return a matrix containing all samples for each observation if \code{summary = FALSE}, else if \code{summary = TRUE} a \code{data.frame} containing a summary of the posterior distribution at each observation (of, respectively, residuals, fitted values, or the spatial trend).
#'
#' \code{plot} returns a \code{ggplot} object that can be customized using the \code{ggplot2} package.
#'
#' \code{as.matrix}, \code{as.data.frame}, \code{as.array} return samples from the joint posterior distribution of parameters in the format corresponding to their names. The \code{pars} argument is used to return samples from only a subset of parameters.
#'
#' @seealso \link[geostan]{stan_glm}, \link[geostan]{stan_esf}, \link[geostan]{stan_icar}, \link[geostan]{stan_car}
#' 
#' @examples 
#' \dontrun{
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
#'                    refresh = 0
#' )
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
#' # extract data.frame of posterior samples
#' S <- as.data.frame(fit, pars = "alpha_re")
#' }
#' 
#' @export
#' @md
#' @method print geostan_fit
#' @rdname geostan_fit
print.geostan_fit <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3, ...) {
  pars <- "intercept"
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
  cat("\nSpatial method (outcome): ", as.character(x$spatial$method), "\n")
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
  if (object$family$family == "auto_gaussian") {
      C <- object$C      
      R <- resid(object, summary = FALSE, detrend = FALSE)
      rho <- as.matrix(object, pars = "car_rho")
      spatial.samples <- t(sapply(1:nrow(rho), function(i) {
          as.numeric( rho[i] * C %*% R[i,] )
      }))
  } else {      
      spatial.samples <- as.matrix(object, pars = par, ...)
  }
  if (summary) {
    return ( post_summary(spatial.samples) )
  } else {
      return( spatial.samples )
  }
}

