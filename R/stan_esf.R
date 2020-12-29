#' Spatial filtering models
#'
#' @export
#' @description Fit a spatial regression model using eigenvector spatial filtering where the spatial filter is 
#' estimated using the regularized horseshoe prior. 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models are specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept specify the grouping variable here using formula synatax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}.
#' @param C Spatial connectivity matrix which will be used to calculate eigenvectors, residual spatial autocorrelation as well as any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms.
#' @param EV A matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial. If provided, still also provide a spatial weights matrix \code{C} for other purposes.  See \link[geostan]{make_EV} and \link[geostan]{shape2mat}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'@param ME To model observational error (i.e. measurement or sampling error) in any or all of the covariates or offset term, provide a named list. Errors are assigned a Gaussian probability distribution and the `true' covariate vector is assigned a Student's t model with optional spatially varying mean. Elements of the list \code{ME} (by name) may include:
#' \describe{
#' 
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one hundred) provide a vector of zeros and ones indicating which columns are bounded. By default the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of the bounded variables. Defaults to \code{bounds = c(0, 100)}.}
#' \item{spatial}{Logical value indicating if the models for covariates should include a spatially varying mean (using an eigenvector spatial filter). Defaults to \code{spatial = FALSE}. If \code{spatial = TRUE} and you do not provide both \code{ME$prior_rhs} and \code{EV} then you must provide a connectivity matrix \code{C}.}
#' \item{prior_rhs}{Optional prior parameters for the regularized horseshoe (RHS) prior used for the ESF data model; only used if \code{ME$spatial = TRUE}. The RHS prior is used for the eigenvector spatial filter (ESF), as in \link[geostan]{stan_esf}. Must be a named list containing vectors \code{slab_df}, \code{slab_scale}, \code{scale_global}, and \code{varname}. The character vector \code{varname} indicates the order of the other parameters (by name).}
#' \item{offset}{if you have an offset term with measurement error, include a vector of standard errors to the list and assign it the name \code{offset}.}
#' }
#' @param nsa Include eigenvectors representing negative spatial autocorrelation? Default \code{nsa = FALSE}. Ignored if \code{EV} is provided.
#' @param threshold Threshold for eigenvector MC value; eigenvectors with values below threshold will be excluded from the candidate set. Default \code{threshold = 0.25}; ignored if \code{EV} is provided. 
#' @param family The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{student_t()} and \code{poisson(link = "log")}, and \code{binomial(link = "logit")}. 
#' @param p0 Number of eigenvector coefficients expected to be far from zero. If this and \code{prior_rhs} are missing, Chun et al.'s (2016) formula will be used to fill this in; see \link[geostan]{exp_pars}.
#' The value of \code{p0} is used to control the prior degree of sparsity in the model.
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_sigma A vector with degrees of freedom, location and scale parameters for the half-Student's t prior on the residual standard deviation \code{sigma}. To use a half-Cauchy prior set degrees of freedom to one; e.g. \code{prior_sigma = c(1, 0, 3)}.
#' @param prior_rhs A named vector with the degrees of freedom and scale of the zero mean Student's t slab for regularizing large eigenvector coefficients, and the prior scale for the global shrinkage parameter; e.g. \code{prior_rhs = c(slab_df = 15, slab_scale = 5, scale_global = 0.5)}. If \code{p0} is provided and prior_rhs is missing, \code{p0} will be used automatically to calculate \code{scale_global}, \code{slab_df} will be set to 15, and \code{slab_scale} will be set to be weakly informative automatically based on the scale of the data.
#' @param prior_nu Set the parameters for the Gamma prior distribution on the degrees of freedom in the likelihood function when using \code{family = student_t}. Defaults to \code{prior_nu = c(alpha = 3, beta = 0.2)}.
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts (\code{alpha_re}). The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 20, location = 0, scale = 5)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model. Parameters from the RHS prios include \code{tau} (the global shrinkage parameter) and \code{lambda} (the local shrinkage parameter).
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.99} and \code{max_treedepth = 15}.
#' @param silent If \code{TRUE}, suppress printed messages including prior specifications and Stan sampling progress (i.e. \code{refresh=0}). Stan's error and warning messages will still print.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details 
#'  The function returns the spatial filter \code{esf}, 
#'  i.e. the linear combination of eigenvectors representing spatial autocorrelation patterns in the outcome variable. 
#'  The entire posterior distribution of the spatial filter
#'  can be obtained with the following code: \code{post_esf <- spatial(fit, summary = FALSE)} where \code{fit} is
#'  the \code{geostan_fit} object returned by a call to \code{stan_esf}. 
#'  
#'  When \code{family = student_t()}, the parameter \code{nu} in the model refers to the degrees of freedom in the Student's t likelihood function for the data.
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{data}{a data frame containing the model data}
#' \item{EV}{A matrix of eigenvectors created with \code{w} and \code{geostan::make_EV}}
#' \item{C}{The spatial weights matrix used to construct EV}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re},  the random effects (varying intercepts) formula if provided, and 
#' \code{data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("esf") and method ("ESF")}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source 
#'
#' Chun, Y., D. A. Griffith, M. Lee and P. Sinha (2016). "Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters." Journal of Geographical Systems, 18(1), 67-85. \link{10.1007/s10109-015-0225-3}
#'
#' Dray, S., P. Legendre & P. R. Peres-Neto (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modeling, 196(3-4), 483-493.
#' 
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. Spatial Statistics. \link{https://doi.org/10.1016/j.spasta.2020.100450}
#'
#' Griffith, Daniel A., and P. R. Peres-Neto (2006). "Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses." Ecology 87(10), 2603-2613.
#' 
#' Griffith, D., and Y. Chun (2014). Spatial autocorrelation and spatial filtering, Handbook of Regional Science. Fischer, MM and Nijkamp, P. eds.
#'
#' Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In Electronic Journal of Statistics, 11(2):5018-5051. \link{https://projecteuclid.org/euclid.ejs/1513306866}
#' 
#' @examples 
#' library(rstan)
#' library(bayesplot)
#' library(ggplot2)
#' library(sf)
#' options(mc.cores = parallel::detectCores()) 
#' 
#' data(ohio)
#' C <- shape2mat(ohio, "B")
#' fit <- stan_esf(gop_growth ~ historic_gop + log(pop_density) + college_educated,
#'                 slx = ~ historic_gop + log(pop_density) + college_educated,
#'                 data = ohio,
#'                 C = C,
#'                 silent = TRUE,
#'                 chains = 1,
#'                 iter = 500,
#'                 family = student_t())
#'
#' # use rstan for sampler diagnostics
#' stan_rhat(fit$stanfit)
#' # see spatial autocorrelation diagnostics
#' spdiag(fit, ohio)
#' 
#' ## summary of estimates
#' plot(fit, pars = c("beta", "intercept", "sigma"))
#' 
#' ## map the spatial filter
#' ohio$esf <- spatial(fit)$mean
#' ggplot(ohio) +
#'       geom_sf(aes(fill = esf)) +
#'       scale_fill_gradient2() +
#'       theme_void()
#' 
#'
stan_esf <- function(formula, slx, re, data, C, EV, ME = NULL,
                     nsa = FALSE, threshold = 0.25,
                     family = gaussian(),
                     p0, prior = NULL, prior_intercept = NULL, prior_sigma = NULL, prior_rhs = NULL, prior_nu = NULL,
                     prior_tau = NULL,
                     centerx = FALSE, scalex = FALSE,
                     prior_only = FALSE,
                     chains = 4, iter = 2e3, refresh = 500, pars = NULL,
                     control = list(adapt_delta = .99, max_treedepth = 15),
                     silent = FALSE,
                     ...) {
  if (missing(formula)) stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (class(family) != "family" | !family$family %in% c("gaussian", "student_t", "poisson", "binomial")) stop ("Must provide a valid family object: gaussian(), student_t(), or poisson().")
  if (missing(C) | missing(data)) stop ("Must provide data and a spatiall connectivity matrix C.")
  if (scalex) centerx <- TRUE
  if (silent) refresh = 0
    ## GLM STUFF -------------  
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")  
  n <- nrow(mod.mat)
  family_int <- family_2_int(family)
  intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0) 
    ## ESF STUFF -------------    
  if (missing(EV)) EV <- make_EV(C, nsa = nsa, threshold = threshold)  
  dev <- ncol(EV)
  if (nrow(EV) != n) stop("nrow(EV) must equal the number of observations.")
    ## GLM STUFF -------------  
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
            if (scalex) Wx <- scale(Wx)                     
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
  ## ESF STUFF -------------
  if (family$family %in% c("poisson", "binomial")) rhs_scale_global = 1
  if (family$family %in% c("gaussian", "student_t")) {
      # if no prior is provided at all for the rhs:
      if (is.null(prior_rhs)) {
          if (missing(p0)) {
              if (missing(C)) stop("To calculate prior_rhs, please provide connectivity matrix C.")
              p0 <- exp_pars(formula = formula, data = tmpdf, C = C)
              if(p0 >= ncol(EV)) p0 <- ncol(EV) - 0.1
              }
       rhs_scale_global <- min(1, p0/(dev-p0) / sqrt(n))
      }
  }
  if (!is.null(prior_rhs)) rhs_scale_global <- as.numeric(prior_rhs["scale_global"])
  ## PARAMETER MODEL STUFF -------------  
  is_student <- family$family == "student_t"
  user_priors <- list(intercept = prior_intercept, beta = prior, sigma = prior_sigma, nu = prior_nu, rhs = prior_rhs, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = user_priors, y = y, x = x, xcentered = centerx,
                        rhs_scale_global = rhs_scale_global, link = family$link, EV = EV, offset = offset)
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
  ## esf data -------------
    dev = dev,
    EV = EV,
    slab_scale = priors$rhs["slab_scale"],
    slab_df = priors$rhs["slab_df"],
    scale_global = priors$rhs["scale_global"],
    prior_only = prior_only
    )
  ## DATA MODEL STUFF -------------  
  me.list <- prep_me_data(ME, x.list$x)
  standata <- c(standata, me.list)
  if (missing(C)) C <- NA
  if (missing(EV)) EV <- NA
  sp.me <- prep_sp_me_data(ME, me.list, C, EV, x.list$x, stan_esf = TRUE, silent = silent)  
  standata <- c(standata, sp.me)
  # -------------  
  # handling multiple possible data types  
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'esf', 'beta_ev', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
  if (is_student) pars <- c(pars, "nu")
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  if (any(me.list$offset_me != 0)) pars <- c(pars, "offset_est")
  priors <- priors[which(names(priors) %in% c(pars, "rhs"))]
    ## PRINT STUFF -------------    
  if (!silent) print_priors(user_priors, priors)
  ## CALL STAN -------------    
   samples <- rstan::sampling(stanmodels$esf, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  if (missing(C)) C <- NA
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x, me.list$x_me_unbounded_idx, me.list$x_me_bounded_idx)  
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$C <- C
  out$EV <- EV  
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME
  out$spatial <- data.frame(par = "esf", method = "ESF")
  class(out) <- append("geostan_fit", class(out))
  return (out)
}


