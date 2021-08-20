#' geostan_fit methods
#'
#' @description Methods for fitted geostan models: extract residuals, fitted values, posterior predictive distribution or spatial component from a spatial regression model; extract samples from the posterior distribution; print regression results; plot posterior distributions.
#'
#' @rdname geostan_fit
#' @param object A fitted model object of class \code{geostan_fit}.
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param summary Logical; should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples from the posterior distribution at each observation is returned.
#' @param pars parameters to include; a character string or vector of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include histograms ("hist"), MCMC traceplots ("trace"), and density plots ("dens"). Diagnostic plots are also available such as Rhat statistics ("rhat"), effective sample size ("ess"), and MCMC autocorrelation ("ac").
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' @param fill fill color for histograms and density plots.
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
#' # extract matrix of samples from posterior distribution of parameters
#' ## alpha_re are the unstructured area random effects
#' S <- as.matrix(fit, pars = "alpha_re")
#'
#' # extract data.frame of posterior samples
#' S <- as.data.frame(fit, pars = "alpha_re")
#' }
#' 
#' @export
#' 
#' @method print geostan_fit
#' @name geostan_fit
#' 
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
  cat("Data models (ME): ")
  if (inherits(x$ME, "list")) {
      if ("se" %in% names(x$ME)) cat(paste(names(x$ME$se), sep = ", "))
      if ("spatial" %in% names(x$ME)) cat("\nData model (prior): CAR (auto Gaussian)") else cat("\nData model (prior): Studen's t")
  } else cat("none")
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
  if (x$spatial$method == "ESF") {
    cat("RHS global shrinkage prior: ", round(x$priors$rhs["scale_global"], 2), "\n")
  }
  print(x$stanfit, pars = pars, digits = digits, probs = probs, ...)
}


#' @export
#' @import graphics
#' @name geostan_fit
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
  rstan::plot(x$stanfit, pars = pars, fill = fill, plotfun = plotfun, ...)
}


#' @export
#' @rdname geostan_fit
#' @method as.matrix geostan_fit
as.matrix.geostan_fit <- function(x, ...) {
  as.matrix(x$stanfit, ...)
}

#' @export
#' @name geostan_fit
#' @method as.data.frame geostan_fit
as.data.frame.geostan_fit <- function(x, ...) {
    as.data.frame(x$stanfit, ...)
}

#' @export
#' @name geostan_fit
#' @method as.array geostan_fit
as.array.geostan_fit <- function(x, ...){
    as.array(x$stanfit, ...)
}

#' @export
#' @method residuals geostan_fit
#' @name geostan_fit
#' 
residuals.geostan_fit <- function(object, summary = TRUE, ...) {
  if (summary) {
    post_summary(object, "residual") 
    } else {
    as.matrix(object$stanfit, pars = "residual")
  }
}

#' @export
#' @import stats
#' @method fitted geostan_fit
#' @name geostan_fit
fitted.geostan_fit <- function(object, summary = TRUE, ...) {
  if (summary) {
    fitted <- post_summary(object, "fitted")
  } else {
   fitted <- as.matrix(object$stanfit, pars = "fitted")
  }
  return(fitted)
}

#' Extract spatial component from a fitted geostan model
#' 
#' @description Extracts the posterior distribution of the spatial component from a fitted geostan model 
#' @seealso \link[geostan]{spatial.geostan_fit}
#' @param object Fitted geostan model
#' @param summary should the posterior distribution be summarized? If \code{FALSE}, returns a matrix of samples; else a \code{data.frame} with summary statistics of the spatial filter at each observation.
#' @param ... additional arguments
#' @export
#' 
spatial <- function(object, summary = TRUE, ...) {
    UseMethod("spatial", object)
}

#' @export
#' @method spatial geostan_fit
#' @name geostan_fit
spatial.geostan_fit <- function(object, summary = TRUE, ...) {
  if (is.na(object$spatial$par)) stop("This model does not have a spatial trend component to extract.")
  par <- as.character(object$spatial$par)
  if (summary) {
    post_summary(object, par)
  } else {
    as.matrix(object$stanfit, pars = par, ...)
  }
}

