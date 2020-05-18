#' geostan_fit methods
#'
#' @description Methods for fitted geostan models: extract residuals, fitted values, posterior predictive distribution or spatial component from a spatial regression model; extract samples from the posterior distribution; print regression results; plot posterior distributions.
#'
#' @rdname geostan_fit
#' @param object A fitted model object of class \code{geostan_fit}.
#' @param x A fitted model object of class \code{geostan_fit}.
#' @param draws Number of samples to return; maximum and default value is the number of MCMC samples stored by the model.
#' @param summary Logical (default \code{summary = TRUE}); should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples from the posterior distribution at each observation is returned.
#' @param pars parameters to include; a character string or vector of parameter names.
#' @param plotfun Argument passed to \code{rstan::plot}. Options include "hist", "trace", and "dens".
#' @param digits number of digits to print
#' @param probs Argument passed to \code{quantile}; which quantiles to calculate and print.
#' @param ... additional arguments.
#' @return return either a matrix containing all samples for each observation, or a \code{data.frame} containing a summary values for each observation.
#' @seealso \link[geostan]{stan_esf} \link[rstan]{stan_plot}
#' @examples 
#' \dontrun{
#' library(geostan)
#' library(sf)
#' data(ohio)
#' fit <- stan_esf(gop_growth ~ historic_gop + log(pop_density), data = ohio, C = shape2mat(ohio), chains = 1, iter = 600)
#'
#' # print and plot results
#' print(fit)
#' plot(fit)
#' plot(fit, pars = c("alpha", "beta")
#' plot(fit, pars = c("b_historic_gop", "b_log(pop_density)"), plotfun = "hist", bins = 30)
#' plot(fit, plotfun = "trace")
#'
#' # extract posterior samples
#' as.data.frame(fit)
#' as.data.frame(fit, pars = "beta")
#' as.data.frame(fit, pars = "b_historic_gop")
#' as.array(fit)
#' as.matrix(fit)
#' 
#' # extract residuals, etc.
#' res <- residuals(fit)
#' ##res <- resid(fit) # alias
#' head(res)
#' plot(ohio$gop_growth, res$mean)
#' mat <- resid(fit, summary = FALSE)
#' dim(mat)
#' pp <- posterior_predict(fit, summary  = FALSE)
#' dim(pp)
#' head(fitted(fit))
#' 
#' # extract and plot the mean estimate of the spatial filter
#' ohio$sf <- spatial(fit, summary = TRUE)$mean
#' plot(ohio[,'sf'])
#' }


#' @export
#' @method print geostan_fit
#' @name geostan_fit
print.geostan_fit <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3, ...) {
  all_pars <- names(x$stanfit)
  pars <- all_pars[which(all_pars %in% c("intercept", "alpha_tau"))]
  if(any(grepl("b_", all_pars))) pars <- c(pars, "beta")
  pars <- c(pars, all_pars[which(all_pars %in% c("sigma", "nu", "p", "rho"))])
  cat("Spatial Regression Results \n")
  cat("Formula: ")
  print(x$formula)
  if (class(x$slx) == "formula") {
   cat("SLX: ")
   print(x$slx)    
  }
  if(!any(is.na(x$re))) {
    cat("Random effects: ")
    print(x$re$formula)
  }
  cat("Spatial method: ", as.character(x$spatial$method), "\n")
  cat("Family: ", x$family$family, "\n")
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
plot.geostan_fit <- function(x, pars, plotfun = "stan_plot", ...) {
  if(missing(pars)) {
    all_pars <- names(x$stanfit)
    pars <- all_pars[which(all_pars %in% c("intercept", "alpha_tau"))]
    if(any(grepl("b_", all_pars))) pars <- c(pars, "beta")
    pars <- c(pars, all_pars[which(all_pars == "sigma")])
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

#' Draw samples from posterior predictive distribution of a fitted geostan model
#' 
#' @description Draws samples from the posterior predictive distribution of a fitted geostan model
#' @param object A fitted geostan model
#' @param draws Number of samples to return; maximum and default value is the number of MCMC samples stored by the model.
#' @param summary Logical (default \code{summary = TRUE}); should the values be summarized with the mean, standard deviation and quantiles (\code{probs = c(.025, .2, .5, .8, .975)}) for each observation? Otherwise a matrix containing samples for each observation is returned.
#' @param ... additional arguments.
#' @export
posterior_predict <- function(object, draws, summary = FALSE, ...) {
    UseMethod("posterior_predict", object)
    }

#' @export
#' @name geostan_fit
#' @method posterior_predict geostan_fit
posterior_predict.geostan_fit <- function(object, draws, summary = FALSE, ...) {
  yrep <- as.matrix(object, pars = "yrep")
  if(!missing(draws)) {
    max_draws <- nrow(yrep)
    idx <- sample(max_draws, size = min(draws, max_draws))
    yrep <- yrep[idx, ]
  }
  if(summary) {
    yrep <- data.frame(
      mean = apply(yrep, 2, mean),
      sd = apply(yrep, 2, sd),
      q.025 = apply(yrep, 2, quantile, probs = .025),
      q.20 = apply(yrep, 2, quantile, probs = .20),
      q.50 = apply(yrep, 2, quantile, probs = .50),
      q.80 = apply(yrep, 2, quantile, probs = .80),
      q.975 = apply(yrep, 2, quantile, probs = .975)
    )
  }
  return(yrep)
}


