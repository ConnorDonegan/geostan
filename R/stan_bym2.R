#' BYM2
#'
#' @export
#' @description Fit the scaled Besag-York-Mollie (BYM2) model introduced by Riebler et al. (2016) with Stan code from Morris et al. (2019). Only fully connected graphs are currenlty supported (i.e. all polygons must have at least one neighbor and there can be no disconnected islands or regions).
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. If an offset term is provided for a Poisson model, it will be transformed to the log scale internally; to add an offset term \code{E} to a formula use \code{y ~ offset(E)}. Binomial models [not yet implemented for \code{stan_bym2}] are specified by setting the left-hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param scaleFactor The scaling factor for the ICAR random effect. Currently INLA is required to calculate this. 
#' @param re If the model includes an additional varying intercept term specify the grouping variable here using formula synatax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param C Spatial connectivity matrix which will be used to construct an edge list, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param family The likelihood function for the outcome variable. Current options are \code{family = poisson(link = "log")}, the default. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. When setting this prior, keep in mind that if \code{centerx = TRUE} (the default), then the intercept is the expected outcome when covariates are at their mean level.
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts (not the exchangeable component of the convolved random effects term, but any additional terms specified by the user). The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 20, location = 0, scale = 20)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{TRUE}. This alters the interpretation of the intercept term, see \code{Details}) below. It also makes setting the prior distribution for the interecept intuitive.
#' @param scalex Should the covariates be scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 5000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details If the \code{centerx = TRUE} (the default), then the intercept is the expected value of the outcome variable when 
#'   all of the covariates are at their mean value.
#' 
#'  The Stan code for the model follows Morris et al. (2019). The \code{INLA} package is currently required to calculate the scaling factor for the spatial random effects. To install it see \code{http://www.r-inla.org/download}. 
#'    
#'  The function returns the ICAR spatial component in the parameter \code{ssre} (spatially structured random effect) and the exchangeable random effects in \code{sure} (spatially unstructured random effect); both parameters are returned after being scaled by \code{sigma_re}. The convolved random effect term (already scaled by \code{sigma_re}) is stored in the parameter \code{convolved_re}. To extract a summary of the posterior distribution for the convolved random effect term from a model, use \code{spatial(fit)}, and for the posterior samples use \code{spatial(fit, summary = FALSE)}.
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the random effects (varying intercepts) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("convolved_re") and method ("BYM2")}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source
#' Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society: Series B (Methodological), 36(2), 192-225.
#'
#' Besag, J., York, J., & Mollié, A. (1991). Bayesian image restoration, with two applications in spatial statistics. Annals of the institute of statistical mathematics, 43(1), 1-20.
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. Spatial and spatio-temporal epidemiology, 31, 100301.
#' 
#' Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical methods in medical research, 25(4), 1145-1165.
#'
stan_bym2 <- function(formula, slx, scaleFactor, re, data, C, family = poisson(),
                     prior = NULL, prior_intercept = NULL,  prior_tau = NULL, 
                centerx = TRUE, scalex = FALSE, chains = 4, iter = 5e3, refresh = 500, pars = NULL,
                control = list(adapt_delta = .9, max_treedepth = 15), ...) {
  if (class(family) != "family" | !family$family %in% c("poisson")) stop ("Must provide a valid family object: poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ offset(E) + x or y ~ 1 for intercept only.")
  if (missing(data) | missing(C)) stop("Must provide data (a data.frame or object coercible to a data.frame) and connectivity matrix C.")
  tmpdf <- as.data.frame(data) 
  nbs <- edges(C)
  n_edges <- nrow(nbs)
  intercept_only <- ifelse(all(dimnames(model.matrix(formula, tmpdf))[[2]] == "(Intercept)"), 1, 0) 
  if (intercept_only) {
    x <- model.matrix(~ 0, data = tmpdf) 
    dx <- 0
    slx <- " "
    scale_params <- list()    
      } else {
    x <- remove_intercept(model.matrix(formula, data = tmpdf))
    if (missing(slx)) {
        slx <- " "
        } else {
           Wx <- SLX(f = slx, DF = tmpdf, SWM = C)
           x <- cbind(Wx, x)
    } 
    dx <- ncol(x)
    x.list <- scale_x(x, center = centerx, scale = scalex)
    x <- x.list$x
    scale_params <- x.list$params    
      }
  ModData <- make_data(formula, tmpdf, x)
  frame <- model.frame(formula, tmpdf)
  y <- model.response(frame)
  n <- nrow(frame)
  if (is.null(model.offset(frame))) {
    offset <- rep(0, times = n)
  } else {
    offset <- model.offset(frame)
    if (family$family == "poisson") offset <- log(offset)
  }
  if(missing(re)) {
    has_re <- n_ids <- id <- 0;
    id_index <- to_index(id, n = nrow(tmpdf))
    re_list <- NA
  } else {
    if (class(re) != "formula") stop("re must be of class formula")
    has_re <- 1
    id <- tmpdf[,paste(re[2])]
    n_ids <- length(unique(id))
    id_index <- to_index(id)
    re_list <- list(formula = re, data = id_index)
  } 
  is_student <- family$family == "student_t"
  priors <- list(intercept = prior_intercept, beta = prior, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = priors, y = y, x = x, xcentered = centerx, link = family$link)
  standata <- list(
    y = y,
    x = x,
    n = n,
    n_edges = n_edges,
    node1 = nbs$node1,
    node2 = nbs$node2,
    dx = dx,
    dim_beta_prior = max(1, dx),
    offset = offset,
    has_re = has_re,
    n_ids = n_ids,
    id = id_index$idx,
    alpha_prior = priors$intercept,
    beta_prior = t(priors$beta), 
    sigma_prior = priors$sigma,
    alpha_tau_prior = priors$alpha_tau,
    phi_scale_prior = 1,
    scaling_factor = scaleFactor,
    is_student = is_student,
    t_nu_prior = priors$nu,    
    has_sigma = family$family %in% c("gaussian", "student_t")
  )
  if (family$family == "binomial") { # not yet implemented
      standata$y <- y[,1]
      standata$N <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'ssre', 'sure', 'sigma_re', 'convolved_re', 'rho', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  priors <- priors[which(names(priors) %in% pars)]
  samples <- rstan::sampling(stanmodels$bym2, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, init_r = 1, ...)
  out <- clean_results(samples, pars, is_student = FALSE, has_re, C, x)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$edges <- nbs
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params  
  out$spatial <- data.frame(par = "convolved_re", method = "BYM2")
  class(out) <- append("geostan_fit", class(out))
  return(out)
}

