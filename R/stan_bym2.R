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
#' @param ME To model measurement error or sampling error in any or all of the covariates, provide a named list containing a dataframe (named \code{ME}) with standard errors for each observation; these will be matched to the variables using column names. If any of the variables in \code{ME} are percentages (ranging from zero to one-hundred---not zero to one!), also include a vector indicating which columns are percentages. For example, if \code{ME} has three columns and the second column is a percentage, include \code{percent = c(0, 1, 0)}. Altogether, \code{ME = list(ME = se.df, percent = c(0, 1, 0))}. This will ensure that the ME models for percentages are properly constrained to the range [0, 100]. Finally, if you have an offset term with measurement error, include a vector of standard errors to the list and assign it the name \code{offset}. The model for offset values will be restricted to allow values only greater than or equal to zero. Note that the \code{ME} model will not work if \code{formula} includes any functions of your variables such as polynomials, splines, or log transformations.
#' @param C Spatial connectivity matrix which will be used to construct an edge list, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param family The likelihood function for the outcome variable. Current options are \code{family = poisson(link = "log")}, the default. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts (not the exchangeable component of the convolved random effects term, but any additional terms specified by the user). The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 20, location = 0, scale = 20)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details 
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
stan_bym2 <- function(formula, slx, scaleFactor, re, data, ME, C, family = poisson(),
                     prior = NULL, prior_intercept = NULL,  prior_tau = NULL, 
                centerx = TRUE, scalex = FALSE, chains = 4, iter = 2e3, refresh = 500, pars = NULL,
                control = list(adapt_delta = .9, max_treedepth = 15), ...) {
  if (class(family) != "family" | !family$family %in% c("poisson")) stop ("Must provide a valid family object: poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ offset(E) + x or y ~ 1 for intercept only.")
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
  is_student = FALSE #// IAR model not compatible with Gaussian/Student's t likelihood
  priors <- list(intercept = prior_intercept, beta = prior, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = priors, y = y, x = x, xcentered = centerx,
                        link = family$link)
  ## DATA MODEL STUFF -------------
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
          if  (!all(names(ME$ME) %in% names(as.data.frame(x)))) stop("All column names in ME$ME must be found in the model matrix (from model.matrix(formula, data)). This error may occur if you've included some kind of data transformation in your model formula, such as a logarithm or polynomial, which is not supported.")
          if (length(ME$percent)) {
              if (length(ME$percent) != ncol(ME$ME)) stop("ME mis-specified: percent must be a vector with one element per column in the ME dataframe.")
              percent <- which(ME$percent == 1)
              not.percent <- which(ME$percent != 1)
          } else {
              percent <- rep(0, times = ncol(ME$ME))
              not.percent <- rep(1, times = ncol(ME$ME))
          }           
                                        # gather any/all variables without ME
          x.df <- as.data.frame(x.list$x)  
          x_obs_idx <- as.array( which( !names(x.df) %in% names(ME$ME) )) 
          x_obs <- as.data.frame(x.df[, x_obs_idx])
          dx_obs <- ncol(x_obs)
                                        # now get all of those with ME
          X.me <- data.frame( x.df[, names(ME$ME)] )
          names(X.me) <- names(ME$ME)
                                        # now X.me needs to be parsed into proportion/non-proportion variables (expressed as percentages)
          x_me_cont <- as.matrix(X.me[,not.percent], nrow = n)
          x_me_prop <- as.matrix(X.me[,percent], nrow = n)
          dx_me_cont <- ncol(x_me_cont)
          dx_me_prop <- ncol(x_me_prop)
                                        # identify which columns in the design matrix correspond to each type of ME variable
          x_me_cont_idx <- as.array( which( names(x.df) %in% names(ME$ME)[not.percent] ))
          x_me_prop_idx <- as.array( which( names(x.df) %in% names(ME$ME)[percent] )) 
          sigma_me_cont <- as.matrix(ME$ME[,not.percent], nrow = n)
          sigma_me_prop <- as.matrix(ME$ME[,percent], nrow = n)
                                        # handle unused parts
          if (!dx_obs) {
              x_obs <- model.matrix(~ 0, tmpdf) 
              x_obs_idx <- a.zero
          }
          if (!dx_me_prop) {
              sigma_me_prop <- x_me_prop <- matrix(0, nrow = n, ncol = 1)
              x_me_prop_idx <- a.zero
          }
          if (!dx_me_cont) {
              sigma_me_cont <- x_me_cont <- matrix(0, nrow = n, ncol = 1)
              x_me_cont_idx <- a.zero
          }
                      # return items in data list ready for Stan: with ME model for covariates
          me.x.list <- list(
          dx_obs = dx_obs,
          dx_me_cont = dx_me_cont,
          dx_me_prop = dx_me_prop,
          x_obs_idx = x_obs_idx,
          x_me_prop_idx = x_me_prop_idx,
          x_me_cont_idx = x_me_cont_idx,
          x_obs = x_obs,
          x_me_prop = x_me_prop,
          x_me_cont = x_me_cont,
          sigma_me_prop = sigma_me_prop,
          sigma_me_cont = sigma_me_cont
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
          dx_me_cont = 0,
          dx_me_prop = 0,
          x_obs_idx = x_obs_idx,
          x_me_prop_idx = a.zero,
          x_me_cont_idx = a.zero,
          x_obs = x_obs,
          x_me_prop = matrix(0, nrow = n, ncol = 1),
          x_me_cont = matrix(0, nrow = n, ncol = 1),
          sigma_me_prop = matrix(0, nrow = n, ncol = 1),
          sigma_me_cont = matrix(0, nrow = n, ncol = 1)
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
          dx_me_cont = 0,
          dx_me_prop = 0,
          x_obs_idx = x_obs_idx,
          x_me_prop_idx = a.zero,
          x_me_cont_idx = a.zero,
          x_obs = x_obs,
          x_me_prop = matrix(0, nrow = n, ncol = 1),
          x_me_cont = matrix(0, nrow = n, ncol = 1),
          sigma_me_prop = matrix(0, nrow = n, ncol = 1),
          sigma_me_cont = matrix(0, nrow = n, ncol = 1),
          offset_me = rep(0, times = n),
          model_offset = 0
      )
  }
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
  ## bym2 data -------------
    n_edges = n_edges,
    node1 = nbs$node1,
    node2 = nbs$node2,
    phi_scale_prior = 1,
    scaling_factor = scaleFactor
    )
  standata <- c(standata, me.list)
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1] ## must overwrite y with vector
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'ssre', 'sure', 'sigma_re', 'convolved_re', 'rho', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  priors <- priors[which(names(priors) %in% pars)]
  samples <- rstan::sampling(stanmodels$bym2, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, init_r = 1, ...)
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$edges <- nbs
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME  
  out$spatial <- data.frame(par = "convolved_re", method = "BYM2")
  class(out) <- append("geostan_fit", class(out))
  return(out)
}

