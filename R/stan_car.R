#' Conditional autoregressive models
#'
#' @export
#' @description Fit a regression model with a conditional auto-regressive (CAR) spatial component. Only binary weighting schemes are currenlty supported.
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized version of the user-provided spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term (or "spatially unstructured random effect") specify the grouping variable here using formula synatax, as in \code{~ ID}. 
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param ME To model observational error (i.e. measurement or sampling error) in any or all of the covariates or offset term, provide a named list. Errors are assigned a Gaussian probability distribution and the `true' covariate vector is assigned a Student's t model with optional spatially varying mean. Elements of the list \code{ME} (by name) may include:
#' \describe{
#' 
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one hundred) provide a vector of zeros and ones indicating which columns are bounded. By default the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of the bounded variables. Defaults to \code{bounds = c(0, 100)}.}
#' \item{spatial}{Logical value indicating if the models for covariates should include a spatially varying mean (using an eigenvector spatial filter). Defaults to \code{spatial = FALSE}. If \code{spatial = TRUE} and you do not provide both \code{ME$prior_rhs} and \code{EV} then you must provide a connectivity matrix \code{C}.}
#' \item{prior_rhs}{Optional prior parameters for the regularized horseshoe (RHS) prior used for the ESF data model; only used if \code{ME$spatial = TRUE}. The RHS prior is used for the eigenvector spatial filter (ESF), as in \link[geostan]{stan_esf}. Must be a named list containing vectors \code{slab_df}, \code{slab_scale}, \code{scale_global}, and \code{varname}. The character vector \code{varname} indicates the order of the other parameters (by name).}
#' \item{offset}{if you have an offset term with measurement error, include a vector of standard errors to the list and assign it the name \code{offset}.}
#' }
#' @param C Binary spatial connectivity matrix. This is used for the CAR prior and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms (it will be row-standardized before calculating \code{slx} terms).
#' @param EV A matrix of eigenvectors from any (transformed) connectivity matrix (presumably spatial). See \link[geostan]{make_EV} and \link[geostan]{shape2mat}. 
#' @param family The likelihood function for the outcome variable. Current options are \code{binomial(link = "logit")} and \code{poisson(link = "log")}. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of varying intercepts \code{re}. The varying intercepts are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 15, location = 0, scale = 5)}.
#' @param prior_phi_precision Hyperprior parameters for the precision (inverse scale) of the CAR component \code{phi}, a numeric vector of length two. The precision parameter for \code{phi} is \code{phi_tau} which is assigned a Gamma prior distribution; defaults to \code{prior_phi_precision = c(2, 2)}.
#' @param centerx Logical value indicating if the covariates should be centered prior to fitting the model.
#' @param scalex Logical value indicating if the covariates be centered and scaled (divided by their standard deviation).
#' @param prior_only Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.
#' @param chains Number of MCMC chains to estimate. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param silent If \code{TRUE}, suppress printed messages including prior specifications and Stan sampling progress (i.e. \code{refresh=0}). Stan's error and warning messages will still print.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details
#' 
#'  The Stan code for the CAR component of the model is from Joseph (2016).
#'   
#'  The CAR model is the prior distribution for the parameter vector \code{phi} and has two associated parameters: \code{phi_alpha} which controls the degree of spatial autocorrelation (and thus the amount of spatial smoothing) and \code{phi_tau} which is a precision (inverse scale) parameter.
#'  
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame, including \code{phi} the spatial term, \code{phi_tau} the scale of the CAR model, and \code{phi_alpha} the parameter that controls the degree of spatial autocorrelation in \code{phi} and has range \code{[0, 1]}.}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors (a sparse representation of C)}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including CAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the random effects (varying intercepts) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("phi") and method ("CAR")}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source
#'
#' Joseph, Max (2016). Exact Sparse CAR Models in Stan. Stan Case Studies, Vol. 3. \url{https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html}
#' 
#' Lawson, Andrew B. (2013). Bayesian Disease Mapping: Hierarchical Modeling in Spatial Epidemiology. CRC Press.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(bayesplot)
#' library(sf)
#' options(mc.cores = parallel::detectCores())
#' data(sentencing)
#'
#' C <- shape2mat(sentencing)
#' log_e <- log(sentencing$expected_sents)
#' fit.car <- stan_car(sents ~ offset(log_e),
#'                     family = poisson(),
#'                     data = sentencing,
#'                     C = C,
#'                     refresh = 0
#' )
#' 
#' # posterior predictive check: predicted distribution should resemble observed distribution
#' yrep <- posterior_predict(fit.car, samples = 75)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep)
#'
#' sp.trend <- spatial(fit.car)$mean
#' st_as_sf(sentencing) %>%
#'   ggplot() +
#'   geom_sf(aes(fill = sp.trend)) +
#'   scale_fill_gradient2()
#'  }
#' 
stan_car <- function(formula, slx, re, data, ME = NULL, C, EV,
                     family = poisson(),
                     prior = NULL, prior_intercept = NULL, prior_tau = NULL, prior_phi_precision = c(2, 2),
                     centerx = FALSE, scalex = FALSE, prior_only = FALSE,
                     chains = 4, iter = 10e3, refresh = 500, pars = NULL,
                     control = list(adapt_delta = .9, max_treedepth = 15),
                     silent = FALSE,
                     ...) {
  if (class(family) != "family" | !family$family %in% c("binomial", "poisson")) stop ("Must provide a valid family object: binomial() or poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(data) | missing(C)) stop("Must provide data (a data.frame or object coercible to a data.frame) and binary connectivity matrix C.")
  if (scalex) centerx = TRUE
  if (silent) refresh = 0
  ## GLM STUFF -------------
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")    
  n <- nrow(mod.mat)
  family_int <- family_2_int(family)
  intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0)
  ## CAR STUFF -------------
  if (any(dim(C) != n)) stop("Dimensions of matrix C must match the number of observations. See ?shape2mat for help creating C.")
  if (!all(C %in% c(0, 1))) stop("C must be a binary connectivity matrix (0s and 1s only). See ?shape2mat for help creating C.")
  nbs <- edges(C)
  n_edges <- nrow(nbs)
  ## GLM STUFF -------------    
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
  user_priors <- list(intercept = prior_intercept, beta = prior, sigma = NULL, nu = NULL, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = user_priors, y = y, x = x, xcentered = centerx,
                        link = family$link, offset = offset)
  ## CAR STUFF -------------
  priors$phi_scale_prior <- prior_phi_precision
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
  ## car data -------------
    C_n = n_edges,
    C_sparse = nbs,
    D_sparse = rowSums(C),
    C = C,
    phi_tau_prior = priors$phi_scale_prior,
  ## draw samples from prior only, ignoring data and likelihood
    prior_only = prior_only
  )
  ## DATA MODEL STUFF -------------  
  me.list <- prep_me_data(ME, x.list$x)
  standata <- c(standata, me.list)
  if (missing(C)) C <- NA
  if (missing(EV)) EV <- NA
  sp.me <- prep_sp_me_data(ME, me.list, C, EV, x.list$x, silent = silent)
  standata <- c(standata, sp.me)
  ## STAN STUFF -------------    
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', "phi", "phi_alpha", "phi_tau", 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  if (any(me.list$offset_me != 0)) pars <- c(pars, "offset_est")
  priors <- priors[which(names(priors) %in% pars)]
  ## PRINT STUFF -------------    
  if (!silent) print_priors(user_priors, priors)
  ## CALL STAN -------------  
  samples <- rstan::sampling(stanmodels$car, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
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
  out$spatial <- data.frame(par = "phi", method = "CAR")
  class(out) <- append("geostan_fit", class(out))
  return (out)
}

