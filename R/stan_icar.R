#' ICAR
#'
#' @export
#' @description Fit a regression model with an intrinsic conditional auto-regressive (ICAR) spatial component. Only fully connected graphs are currenlty supported (i.e. all polygons must have at least one neighbor and there can be no disconnected islands or regions).
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. If an offset term is provided for a Poisson model, it will be transformed to the log scale internally; to add an offset term \code{E} to a formula use \code{y ~ offset(E)}. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
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
#' \item{offset}{if you have an offset term with measurement error, include a vector of standard errors to the list and assign it the name \code{offset}. The model for offset values will be restricted to allow values only greater than or equal to zero.}
#' }
#' @param C Spatial connectivity matrix which will be used to construct an edge list, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param family The likelihood function for the outcome variable. Current options are \code{binomial(link = "logit")} and \code{poisson(link = "log")}. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts. The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 15, location = 0, scale = 5)}.
#' @param prior_phi Prior for the scale of the spatial ICAR component \code{phi}. \code{phi} is scaled by the parameter \code{phi_scale} which is given a positively-constrained (half-) Gaussian prior distribution with its location parameter at zero and scale parameter set to \code{prior_phi}. This defaults to \code{prior_phi = 1}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
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
#'  When \code{family = student_t()}, the parameter \code{nu} in the model refers to the degrees of freedom in the Student's t likelihood function for the data.
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
#' fit.icar <- stan_icar(sents ~ offset(expected_sents), family = poisson(),
#'                      data = sentencing@data, C = C,
#'                     chains = 1, iter = 500)
#'
#' # add exchangeable random effects (for a BYM model). See stan_bym2 for a better model.
#' fit.bym <- stan_icar(sents ~ offset(expected_sents), re = ~ name,
#'                      family = poisson(),
#'                      data = sentencing@data, C = C,
#'                      chains = 1, iter = 500)
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
stan_icar <- function(formula, slx, re, data, ME, C, family = poisson(),
                      prior = NULL, prior_intercept = NULL, prior_tau = NULL, prior_phi = 1,
                centerx = FALSE, scalex = FALSE, chains = 4, iter = 2e3, refresh = 500, pars = NULL,
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
  if (intercept_only) {
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
            dwx <- ncol(Wx)
            dw_nonzero <- sum(W!=0)
            wx_idx <- as.array( which(paste0("w.", dimnames(x)[[2]]) %in% dimnames(Wx)[[2]]), dim = dwx )
            x <- cbind(Wx, x)
    }
    dbeta_prior <- ncol(x) ## dimensions of beta prior; x includes slx, if any; x.list$x is the processed model matrix without slx terms.
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
  # some defaults
  dx_me_unbounded <- 0
  dx_me_bounded <- 0
  x_me_bounded_idx = a.zero
  x_me_unbounded_idx = a.zero
  bounds <- c(0, 100)
  if (!missing(ME)) {
      if (!inherits(ME, "list")) stop("ME must be a list .")
                # ME model for offset
      if (length(ME$offset)) {
          offset_me <- ME$offset
          model_offset <- 1
      } else {
          offset_me <- rep(0, times = n)
          model_offset <- 0
      }
            # return items in data list ready for Stan: with ME model for offset
      me.list <- list(offset_me = offset_me,
                      model_offset = model_offset
                      )
      if (!is.null(ME$ME)) {
          if (!inherits(ME$ME, "data.frame")) stop("ME must be a list in which the element named ME is of class data.frame, containing standard errors for the observations.")
          if  (!all(names(ME$ME) %in% names(as.data.frame(x)))) stop("All column names in ME$ME must be found in the model matrix (from model.matrix(formula, data)). This error may occur if you've included some kind of data transformation in your model formula, such as a logarithm or polynomial, which is not supported for variables with sampling/measurement error.")
          if (length(ME$bounded)) {
              if (length(ME$bounded) != ncol(ME$ME)) stop("ME mis-specified: bounded must be a vector with one element per column in the ME dataframe.")
              bounded <- which(ME$bounded == 1)
              not.bounded <- which(ME$bounded != 1)
              if (length(ME$bounds)) {
                  if(length(ME$bounds) != 2 | !inherits(ME$bounds, "numeric")) stop("ME$bounds must be numeric vector of length 2.")
                  bounds <- ME$bounds
              }
          } else {
              bounded <- integer(0) #rep(0, times = ncol(ME$ME))
              not.bounded <- 1:ncol(ME$ME) #rep(1, times = ncol(ME$ME))
          }           
                                        # gather any/all variables without ME
          x.df <- as.data.frame(x.list$x)  
          x_obs_idx <- as.array( which( !names(x.df) %in% names(ME$ME) )) 
          x_obs <- as.data.frame(x.df[, x_obs_idx])
          dx_obs <- ncol(x_obs)
                                        # now get all of those with ME
          ## X.me <- data.frame( x.df[, names(ME$ME)] )
          ## names(X.me) <- names(ME$ME)
          
                                        # now X.me needs to be parsed into bounded/non-bounded variables and ordered as x
          nm_me_unbounded <- names(ME$ME)[not.bounded]
          x_me_unbounded <- data.frame( x.df[, nm_me_unbounded] )
          names(x_me_unbounded) <- nm_me_unbounded
          x_me_unbounded_order <- na.omit( match(names(x.df), names(x_me_unbounded)) )
          x_me_unbounded <- data.frame(x_me_unbounded[, x_me_unbounded_order])
          dx_me_unbounded <- ncol(x_me_unbounded)
          sigma_me_unbounded <- as.matrix(ME$ME[,not.bounded], nrow = n)
          sigma_me_unbounded <- as.matrix(sigma_me_unbounded[, x_me_unbounded_order], nrow = n)
          x_me_unbounded_idx <- as.array( which( names(x.df) %in% nm_me_unbounded ))
          
          nm_me_bounded <- names(ME$ME)[bounded]
          x_me_bounded <- data.frame( x.df[, nm_me_bounded] )
          names(x_me_bounded) <- nm_me_bounded
          x_me_bounded_order <- na.omit( match(names(x.df), names(x_me_bounded)) )
          x_me_bounded <- as.matrix(x_me_bounded[, x_me_bounded_order], nrow = n)
          dx_me_bounded <- ncol(x_me_bounded)
          sigma_me_bounded <- as.matrix(ME$ME[,bounded], nrow = n)
          sigma_me_bounded <- as.matrix(sigma_me_bounded[, x_me_bounded_order], nrow = n)
          x_me_bounded_idx <- as.array( which( names(x.df) %in% nm_me_bounded ))

                                        # handle unused parts
          if (!dx_obs) {
              x_obs <- model.matrix(~ 0, tmpdf) 
              x_obs_idx <- a.zero
          }
          if (!dx_me_bounded) {
              sigma_me_bounded <- x_me_bounded <- matrix(0, nrow = n, ncol = 1)
              x_me_bounded_idx <- a.zero
          }
          if (!dx_me_unbounded) {
              sigma_me_unbounded <- x_me_unbounded <- matrix(0, nrow = n, ncol = 1)
              x_me_unbounded_idx <- a.zero
          }
                      # return items in data list ready for Stan: with ME model for covariates
          me.x.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,
          x_obs = x_obs,
          x_me_bounded = x_me_bounded,
          x_me_unbounded = x_me_unbounded,
          sigma_me_bounded = sigma_me_bounded,
          sigma_me_unbounded = sigma_me_unbounded
      )
          me.list <- c(me.list, me.x.list)
      } else {
      # in this case there is just an offset me model but not ME.X
      x_obs <- x.list$x
      dx_obs <- ncol(x_obs)
          if (dx_obs) {
              x_obs_idx <- as.array(1:dx_obs, dim = dx_obs)
          } else {
              x_obs_idx <- a.zero
          }
      me.x.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,          
          x_obs = x_obs,
          x_me_bounded = matrix(0, nrow = n, ncol = 1),
          x_me_unbounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_bounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_unbounded = matrix(0, nrow = n, ncol = 1)
      )
     me.list <- c(me.list, me.x.list)
      }   
  }
      # return items in data list ready for Stan: no ME model at all
  if (missing(ME)) {
      x_obs <- x.list$x
      dx_obs <- ncol(x_obs)
      if (dx_obs) {
          x_obs_idx <- as.array(1:dx_obs, dim = dx_obs)
      } else {
          x_obs_idx <- a.zero
      }
      me.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,          
          x_obs = x_obs,
          x_me_bounded = matrix(0, nrow = n, ncol = 1),
          x_me_unbounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_bounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_unbounded = matrix(0, nrow = n, ncol = 1),
          offset_me = rep(0, times = n),
          model_offset = 0
      )
  }
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
    phi_scale_prior = priors$phi_scale_prior
    )
  standata <- c(standata, me.list)
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  ## STAN STUFF -------------    
  pars <- c(pars, 'intercept', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (dx_me_bounded) pars <- c(pars, "x_true_bounded")
  priors <- priors[which(names(priors) %in% pars)]
  ## CALL STAN -------------  
  samples <- rstan::sampling(stanmodels$icar, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, init_r = 1, ...)
  if (missing(C)) C <- NA
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x, x_me_unbounded_idx, x_me_bounded_idx)
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

