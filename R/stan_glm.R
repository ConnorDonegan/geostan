#' GLM
#'
#' @export
#' @description Fit a Poisson generalized linear model with (optional) exchangeable random effects. 
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. If an offset term is provided for a Poisson model, it will be transformed to the log scale (and it will be ignored if its not a Poisson model); to add an offset term \code{E} to a formula use \code{y ~ offset(E)}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally); you must also provide \code{C} when including \code{slx}.
#'  Specified covariates will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term (or "spatially unstructured random effect") specify the grouping variable here using formula synatax, as in \code{~ ID}.  The resulting random effects parameter returned is named \code{alpha_re}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param C Optional spatial connectivity matrix which will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param family The likelihood function for the outcome variable. Current options are \code{family = poisson(link = "log")}. 
#' @param prior A \code{data.frame} or \code{matrix} with Student's t prior parameters for the coefficients. Provide three columns---degrees of freedom, location and scale---and a row for each variable in their order of appearance in the model formula. For now, if you want a Gaussian prior use very large degrees of freedom. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with degrees of freedom, location and scale parameters for a Student's t prior on the intercept; e.g. \code{prior_intercept = c(15, 0, 10)}.
#' @param prior_sigma A vector with degrees of freedom, location and scale parameters for the half-Student's t prior on sigma^2; e.g. \code{prior_sigma = c(5, 0, 10)}. To use a half-Cauchy prior set degrees of freedom to one. 
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts. The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 20, location = 0, scale = 20)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{TRUE} for computational efficiency. This alters the interpretation of the intercept term! See \code{Details}) below.
#' @param scalex Should the covariates be scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 5000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details If the \code{centerx = TRUE} (the default), then the intercept is the expected value of the outcome variable when 
#'   all of the covariates are at their mean value. This often has interpretive value in itself
#'   though it is the  default here for computational reasons.
#'  
#'  When \code{family = student_t()}, the parameter \code{nu} in the model refers to the degrees of freedom in the Student's t likelihood function for the data.
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), residual spatial autocorrelation (Moran coefficient of the residuals), and root mean squared error. Residuals are relative to the mean fitted value for each observation.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the random effects (varying intercepts) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("phi") and method ("ICAR")}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @examples
#' library(ggplot2)
#' library(sf)
#' data(sentencing)
#'
#' # using a small number of iterations and a single chain only for compilation speed
#' # add C to calculate residual SA internally
#' C <- shape2mat(sentencing)
#' fit.pois <- stan_glm(sents ~ offset(expected_sents),
#'                      re = ~ name,
#'                      family = poisson(),
#'                      data = sentencing@data, C = C,
#'                     chains = 1, iter = 2e3)
#'
#' # diagnostics plot: Rhat values should all by very near 1
#' library(rstan)
#' rstan::stan_rhat(fit.pois$stanfit)
#'  # see effective sample size for all parameters and generated quantities
#'  # (including residuals, predicted values, etc.)
#' rstan::stan_ess(fit.pois$stanfit)
#' # or for a particular parameter
#' rstan::stan_ess(fit.pois$stanfit, "phi")

# posterior predictive check: predicted distribution should resemble observed distribution
#' library(bayesplot)
#' yrep <- posterior_predict(fit.pois)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep[1:100,])

#' # map standardized incidence (sentencing) ratios
#' E <- sentencing$expected_sents
#' sentencing@data$ssr <- fitted(fit.pois)$mean / E
#' ggplot(data=st_as_sf(sentencing),
#'       aes(fill = ssr)) +
#'  geom_sf() +
#'  scale_fill_gradient2(midpoint = 1) +
#'  theme_bw() +
#'  ggtitle("Standardized state prison sentencing ratios, 1905-1910")
#'
stan_glm <- function(formula, slx, re, data, C, family = gaussian(),
                      prior = NULL, prior_intercept = NULL, prior_sigma = NULL, prior_nu = NULL,
                      prior_tau = NULL,
                centerx = TRUE, scalex = FALSE, chains = 4, iter = 5e3, refresh = 500, pars = NULL,
                control = list(adapt_delta = .9, max_treedepth = 15), ...) {
  if (class(family) != "family" | !family$family %in% c("gaussian", "student_t", "poisson")) stop ("Must provide a valid family object: poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(data)) stop("Must provide data (a data.frame or object coercible to a data.frame).")
  tmpdf <- as.data.frame(data)
  intercept_only <- ifelse(all(dimnames(model.matrix(formula, tmpdf))[[2]] == "(Intercept)"), 1, 0) 
  if (intercept_only) {
    x <- model.matrix(~ 0, data = tmpdf) 
    dx <- 0
    slx <- " "
      } else {
    xraw <- model.matrix(formula, data = tmpdf)
    xraw <- remove_intercept(xraw)
    x <- scale_x(xraw, center = centerx, scale = scalex)
    if (missing(slx)) {
        slx <- " "
        } else {
           Wx <- SLX(f = slx, DF = tmpdf, SWM = C, cx = centerx, sx = scalex)
           x <- cbind(Wx, x)
    } 
    dx <- ncol(x)
      }
  ModData <- make_data(formula, tmpdf, x)
  frame <- model.frame(formula, tmpdf)
  y <- model.response(frame)
  n <- length(y)
  if (is.null(model.offset(frame))) {
    offset <- rep(0, times = n)
  } else {
    offset <- log(model.offset(frame))
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
    re_list <- list(formula = re, Data = id_index)
  } 
  is_student <- family$family == "student_t"
  priors <- list(intercept = prior_intercept, beta = prior, sigma = prior_sigma, nu = prior_nu, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = priors, y = y, x = x, xcentered = centerx,
                        link = family$link)
    standata <- list(
    y = y,
    x = x,
    n = n,
    n_edges = n_edges,
    node1 = nbs$node1,
    node2 = nbs$node2,
    dx = dx,
    dim_beta_prior = max(1, dx),
    log_E = offset,
    has_re = has_re,
    n_ids = n_ids,
    id = id_index$idx,
    alpha_prior = priors$intercept,
    beta_prior = t(priors$beta), 
    sigma_prior = priors$sigma,
    alpha_tau_prior = priors$alpha_tau,
    t_nu_prior = priors$nu,
    is_student = is_student
    )
  pars <- c(pars, 'intercept', 'phi', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
  if (is_student) pars <- c(pars, "nu")
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  priors <- priors[which(names(priors) %in% pars)]
  ## if (family$family %in% c("gaussian", "student_t")) {
  ##   if (intercept_only) {
  ##     samples <- rstan::sampling(stanmodels$glm_io, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  ##    } else {
  ##      samples <- rstan::sampling(stanmodels$glm_continuous, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  ##    }
  ##   }
  if (family$family == "poisson") {
      samples <- rstan::sampling(stanmodels$glm_poisson, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, init_r = 1, ...)
    }
  if (missing(C)) C <- NA
  out <- clean_results(samples, pars, is_student, has_re, C, x)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$edges <- nbs
  out$re <- re_list
  out$priors <- priors
  out$spatial <- data.frame(par = "alpha_re", method = "Exchangeable")
  class(out) <- append("geostan_fit", class(out))
  return(out)
}

