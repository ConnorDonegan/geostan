#' Intrinsic autoregressive models
#'
#' @export
#' @description Fit a regression model with an intrinsic conditional auto-regressive (ICAR) spatial component. Only fully connected graphs are currenlty supported (i.e. all polygons must have at least one neighbor and there can be no disconnected islands or regions).
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term (or "spatially unstructured random effect") specify the grouping variable here using formula synatax, as in \code{~ ID}. If this is specified at the observational unit level then it becomes the original Besag-York-Mollie model. In that case this random effects term and the ICAR component are not separately identifiable. The resulting random effects parameter returned is named \code{alpha_re}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param ME  To model observational error (i.e. measurement or sampling error) in any or all of the covariates or offset term, provide a named list. Observational errors are assigned a Gaussian probability distribution. Elements of the list \code{ME} (by name) may include:
#' \describe{
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one-hundred) provide a vector of zeros and ones indicating which columns are bounded. By default the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of the bounded variables. Defaults to \code{bounds = c(0, 100)}.}
#' \item{offset}{if you have an offset term with measurement error, include a vector of standard errors to the list and assign it the name \code{offset}.}
#' }
#' @param C Spatial connectivity matrix which will be used to construct an edge list, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param family The likelihood function for the outcome variable. Current options are \code{binomial(link = "logit")} and \code{poisson(link = "log")}. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts. The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 15, location = 0, scale = 5)}.
#' @param prior_phi Prior for the scale of the spatial ICAR component \code{phi}. \code{phi} is scaled by the parameter \code{phi_scale} which is given a positively-constrained (half-) Gaussian prior distribution with its location parameter at zero and scale parameter set to \code{prior_phi}. This defaults to \code{prior_phi = 1}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details
#'  The Stan code for the ICAR component of the model follows Morris et al. (2019).
#'    
#'  The function returns the ICAR component in the parameter \code{phi}.
#'  The entire posterior distribution of \code{phi}
#'  can be obtained with the following code: \code{post_phi <- spatial(fit, summary = FALSE)} where \code{fit} is
#'  the \code{geostan_fit} object returned by a call to \code{stan_icar}. 
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
#' \item{spatial}{A data frame with the name of the spatial component parameter ("phi") and method ("ICAR")}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source
#'
#' Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society: Series B (Methodological), 36(2), 192-225.
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. Spatial and spatio-temporal epidemiology, 31, 100301.
#'
#' @examples
#' library(ggplot2)
#' library(sf)
#' library(rstan)
#' options(mc.cores = parallel::detectCores())
#' data(sentencing)
#'
#' # using a small number of iterations and a single chain only for compilation speed
#' C <- shape2mat(sentencing)
#' log_e <- log(sentencing$expected_sents)
#' fit.icar <- stan_icar(sents ~ offset(log_e),
#'                      family = poisson(),
#'                      data = sentencing,
#'                      C = C,
#'                      cores = 1,  # cores = 4,
#'                      chains = 1, # chains = 4,
#'                     iter = 200)  # iter = 2e3
#'
#' # add exchangeable random effects (for a BYM model). See stan_bym2 for a better model.
#' fit.bym <- stan_icar(sents ~ offset(log_e),
#'                      re = ~ name,
#'                      family = poisson(),
#'                      data = sentencing,
#'                      C = C,
#'                      cores = 1,  # cores = 4,   
#'                      chains = 1, # chains = 4,
#'                      iter = 200) # iter = 2e3
#'
#' # view WAIC 
#' waic(fit.icar)
#' waic(fit.bym)
#'
#'# diagnostics plot: Rhat values should all by very near 1
#' library(rstan)
#' rstan::stan_rhat(fit.icar$stanfit)
#'  # see effective sample size for all parameters and generated quantities
#'  # (including residuals, predicted values, etc.)
#' rstan::stan_ess(fit.icar$stanfit)
#' # or for a particular parameter
#' rstan::stan_ess(fit.icar$stanfit, "phi")

# posterior predictive check: predicted distribution should resemble observed distribution
#' library(bayesplot)
#' yrep <- posterior_predict(fit.icar, samples = 100)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep)

#' # map standardized incidence (sentencing) ratios
#' E <- sentencing$expected_sents
#' sentencing@data$ssr <- fitted(fit.bym)$mean / E
#' ggplot(data=st_as_sf(sentencing),
#'       aes(fill = ssr)) +
#'  geom_sf() +
#'  scale_fill_gradient2(midpoint = 1) +
#'  theme_bw() +
#'  ggtitle("Standardized state prison sentencing ratios, 1905-1910")
#'
stan_icar <- function(formula, slx, re, data, ME = NULL, C, family = poisson(),
                      prior = NULL, prior_intercept = NULL, prior_tau = NULL, prior_phi = 1,
                      centerx = FALSE, scalex = FALSE, prior_only = FALSE,
                      chains = 4, iter = 2e3, refresh = 500, pars = NULL,
                control = list(adapt_delta = .9, max_treedepth = 15), ...) {
  if (class(family) != "family" | !family$family %in% c("binomial", "poisson")) stop ("Must provide a valid family object: binomial() or poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(data) | missing(C)) stop("Must provide data (a data.frame or object coercible to a data.frame) and connectivity matrix C.")
  if (scalex) centerx = TRUE
  ## IAR STUFF -------------
  nbs <- edges(C)
  n_edges <- nrow(nbs)  
  ## GLM STUFF -------------
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  n <- nrow(mod.mat)
  family_int <- family_2_int(family)
  intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0) 
  if (intercept_only) { # x includes slx, if any, for prior specifictions; x.list$x is the processed model matrix without slx terms.
    x <- model.matrix(~ 0, data = tmpdf) 
    dbeta_prior <- 0
    slx <- " "
    scale_params <- list()
    x.list <- list(x = x)
    W <- matrix(0, nrow = 1, ncol = 1)
    dwx <- 0
    dw_nonzero <- 0
    wx_idx <- a.zero
      } else {
    xraw <- model.matrix(formula, data = tmpdf)
    xraw <- remove_intercept(xraw)
    x.list <- scale_x(xraw, center = centerx, scale = scalex)
    x <- x.list$x
    scale_params <- x.list$params
    if (missing(slx)) {
        slx <- " "
        W <- matrix(0, ncol = 1, nrow = 1)
        dwx = 0
        wx_idx = a.zero
        dw_nonzero <- 0
    } else {
            W <- C / rowSums(C)
            Wx <- SLX(f = slx, DF = tmpdf, SWM = W)
            if (scalex) Wx <- scale(Wx)            
            dwx <- ncol(Wx)
            dw_nonzero <- sum(W!=0)
            wx_idx <- as.array( which(paste0("w.", dimnames(x)[[2]]) %in% dimnames(Wx)[[2]]), dim = dwx )
            x <- cbind(Wx, x)
    }
    dbeta_prior <- ncol(x) ## dimensions of beta prior; 
      }
  ModData <- make_data(formula, tmpdf, x)
  frame <- model.frame(formula, tmpdf)
  y <- y_int <- model.response(frame)
  if (family_int %in% c(1,2)) y_int <- rep(0, length(y))
  if (is.null(model.offset(frame))) {
    offset <- rep(0, times = n)
  } else {
      offset <- model.offset(frame)
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
  ## PARAMETER MODEL STUFF -------------  
  is_student <- family$family == "student_t" # always false
  priors <- list(intercept = prior_intercept, beta = prior, sigma = NULL, nu = NULL, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = priors, y = y, x = x, xcentered = centerx,
                        link = family$link)
  ## IAR STUFF -------------
  priors$phi_scale_prior <- prior_phi
  ## DATA MODEL STUFF -------------
  me.list <- prep_me_data(ME, family, x.list$x)
  ## MIXED STUFF -------------  
  standata <- list(
  ## glm data -------------
    y = y,
    y_int = y_int,
    trials = rep(0, length(y)),
    n = n,
    dbeta_prior = dbeta_prior,
    offset_obs = offset,
    has_re = has_re,
    n_ids = n_ids,
    id = id_index$idx,
    alpha_prior = priors$intercept,
    beta_prior = t(priors$beta), 
    sigma_prior = priors$sigma,
    alpha_tau_prior = priors$alpha_tau,
    t_nu_prior = priors$nu,
    family = family_int,
  ## slx data -------------    
    W = W,
    dwx = dwx,
    wx_idx = wx_idx,
  ## iar data -------------
    n_edges = n_edges,
    node1 = nbs$node1,
    node2 = nbs$node2,
    phi_scale_prior = priors$phi_scale_prior,
    prior_only = prior_only
    )
  standata <- c(standata, me.list)
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  ## STAN STUFF -------------    
  pars <- c(pars, 'intercept', "phi", 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  if (any(me.list$offset_me != 0)) pars <- c(pars, "offset_est")
  priors <- priors[which(names(priors) %in% pars)]
  ## CALL STAN -------------  
  samples <- rstan::sampling(stanmodels$icar, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  if (missing(C)) C <- NA
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x, me.list$x_me_unbounded_idx, me.list$x_me_bounded_idx)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$edges <- nbs  
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME
  out$spatial <- data.frame(par = "phi", method = "ICAR")
  class(out) <- append("geostan_fit", class(out))
  return (out)
}

