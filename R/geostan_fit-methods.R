#' geostan_fit methods
#'
#' @description Methods for fitted geostan models: extract residuals, fitted values, posterior predictive distribution or spatial component from a spatial regression model; extract samples from the posterior distribution; print regression results; plot posterior distributions.
#'
#' @rdname geostan_fit
#' @param object A fitted model object of class \code{geostan_fit}.
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param summary Logical; should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples from the posterior distribution at each observation is returned.
#' @param pars parameters to include; a character string or vector of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include "hist", "trace", and "dens".
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' @param ... additional arguments.
#' @return
#'
#' Methods \code{residuals}, \code{fitted}, \code{spatial} return a matrix containing all samples for each observation if \code{summary = FALSE}, else if \code{summary = TRUE} a \code{data.frame} containing a summary of the posterior distribution at each observation (of, respectively, residuals, fitted values, or the spatial trend).
#'
#' \code{plot} returns a \code{ggplot} object that can be customized using the \code{ggplot2} package.
#'
#' \code{as.matrix}, \code{as.data.frame}, \code{as.array} return samples from the (joint) posterior distribution of parameters in the format corresponding to their names. The \code{pars} argument is used to return samples from only a subset of parameters.
#'
#' @seealso \link[geostan]{stan_glm}, \link[geostan]{stan_esf}, \link[geostan]{stan_icar}, \link[geostan]{stan_car}
#' 
#' @examples 
#' \dontrun{
#' library(ggplot2)
#' library(sf)
#' data(ohio)
#' fit <- stan_esf(gop_growth ~ historic_gop + log(pop_density),
#'                data = ohio,
#'                scalex = TRUE,
#'                C = shape2mat(ohio),
#'                refresh = 0 # less printing
#' )
#'
#' # print and plot results
#' print(fit)
#' plot(fit)
#'
#' beta.samples <- as.matrix(fit, pars = "beta")
#' beta.mean <- apply(beta.samples, 2, mean)
#' print(beta.mean)
#' 
#' # extract residuals and fitted values
#' res <- residuals(fit)
#' res <- resid(fit) 
#' head(res)
#' f <- fitted(fit)
#' plot(f$mean, res$mean, xlab = "fitted values", ylab = "residuals")
#' mat <- resid(fit, summary = FALSE)
#' dim(mat)
#' 
#' # extract and plot the posterior mean of the spatial filter
#' sf <- spatial(fit, summary = TRUE)
#' head(sf)
#' ohio$sf <- sf$mean
#' ggplot(ohio) +
#'   geom_sf(aes(fill = sf)) +
#'   scale_fill_gradient2()  
#' }
#' 
#' @export
#' 
#' @method print geostan_fit
#' @name geostan_fit
#' 
print.geostan_fit <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3, ...) {
  pars <- "intercept"
  cat("Spatial Regression Results \n")
  cat("Formula: ")
  print(x$formula)
  if (class(x$slx) == "formula") {
   cat("SLX: ")
   print(x$slx)
   pars <- c(pars, "gamma")
  }
  x.pars <- c("beta", "nu", "sigma")
  if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])   
  if(!any(is.na(x$re))) {
    cat("Partial pooling (varying intercept): ")
    print(x$re$formula)
    pars <- c(pars, "alpha_tau")
  }
  cat("Data models (ME): ")
  if (inherits(x$ME, "list")) {
      if ("se" %in% names(x$ME)) cat(paste(names(x$ME$se), sep = ", "))
      if ("spatial" %in% names(x$ME)) cat("\nPrior data model: Student's t with spatially varying mean (ESF)") else cat("\nPrior data model: Studen's t")
  } else cat("none")
  cat("\nSpatial method (outcome): ", as.character(x$spatial$method), "\n")
  if (x$spatial$method == "CAR") pars <- c(pars, "phi_alpha", "phi_tau")
  if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
  if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
  if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")  
  cat("Likelihood function: ", x$family$family, "\n")
  cat("Link function: ", x$family$link, "\n")
  cat("Residual Moran Coefficient: ", x$diagnostic[grep("Residual_MC", attributes(x$diagnostic)$names)], "\n")
  cat("WAIC: ", x$diagnostic[grep("WAIC", attributes(x$diagnostic)$names)], "\n")
  cat("Observations: ", nrow(x$data), "\n")
  if(x$spatial[1] == "esf") {
    cat("RHS global shrinkage prior: ", round(x$priors$rhs["scale_global"], 2), "\n")
  }
  print(x$stanfit, pars = pars, digits = digits, probs = probs, ...)
}


#' @export
#' @import graphics
#' @name geostan_fit
#' @method plot geostan_fit
plot.geostan_fit <- function(x, pars, plotfun = "dens", ...) {
  if(missing(pars)) {
      pars <- "intercept"
      x.pars <- c("beta", "nu", "sigma")
      if (any(x.pars %in% names(x$priors))) pars <- c(pars, names(x$priors)[grep(paste0(x.pars, collapse="|"), names(x$priors))])
      if (x$spatial$method == "CAR") pars <- c(pars, "phi_alpha", "phi_tau")
      if (x$spatial$method == "BYM2") pars <- c(pars, "rho", "spatial_scale")
      if (x$spatial$method == "BYM") pars <- c(pars, "spatial_scale", "theta_scale")
      if (x$spatial$method == "ICAR") pars <- c(pars, "spatial_scale")
  }
  rstan::plot(x$stanfit, pars = pars, plotfun = plotfun, ...)
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
  par <- as.character(object$spatial$par)
  if (summary) {
    post_summary(object, par)
  } else {
    as.matrix(object$stanfit, pars = par, ...)
  }
}

