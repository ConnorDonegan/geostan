#' Regression with spatial filtering
#'
#' @export
#' @description Fit a spatial regression model using eigenvector spatial filtering where the spatial filter is 
#' estimated using the regularized horseshoe prior. 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. If an offset term is provide for a Poisson model, \code{log(offset)} is passed to the Stan model. (To add an offset term \code{E} to your model formula use \code{y ~ offset(E)}).
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept specify the grouping variable here using formula synatax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}.
#' @param shape A simple features (\code{sf}) object or \code{SpatialPolygonsDataFrame}. A binary neighbors matrix will be constructed by queen contiguity condition and passed to \code{make_EV} to obtain eigenvectors for the model. Alternatively, supply \code{EV} and \code{data}.
#' @param EV Matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial. If provided, also provide spatial weights matrix \code{C}.  See \link[geostan]{make_EV} and \link[geostan]{shape2mat}. This \code{EV} argument is ignored if \code{shape} is provided.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data. Not used if \code{shape} is provided.
#' @param C Optional spatial connectivity matrix which, if provided, will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms in the event that \code{shape} has not been provided; it will be row-standardized before calculating \code{slx} terms.
#' @param nsa Include eigenvectors representing negative spatial autocorrelation? Default \code{nsa = FALSE}
#' @param threshold Threshold for eigenvector MC value; eigenvectors with values below threshold will be excluded from the candidate set. Default \code{threshold = .2}
#' @param family The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{family = student_t()} and \code{family = poisson(link = "log")}. 
#' @param p0 Number of eigenvector coefficients expected to be far from zero. If missing, Chun et al.'s (2016) formula will be used to fill this in; see \link[geostan]{exp_pars}.
#' The value of \code{p0} is used to control the prior degree of sparsity in the model.
#' @param prior A \code{data.frame} or \code{matrix} with Student's t prior parameters for the coefficients. Provide three columns---degrees of freedom, location and scale---and a row for each variable in their order of appearance in the model formula. For now, if you want a Gaussian prior use very large degrees of freedom. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with degrees of freedom, location and scale parameters for a Student's t prior on the intercept; e.g. \code{prior_intercept = c(15, 0, 10)}.
#' @param prior_sigma A vector with degrees of freedom, location and scale parameters for the half-Student's t prior on sigma^2. Use a half-Cauchy prior by setting degrees of freedom to one; e.g. \code{prior_sigma = c(5, 0, 10)}.
#' @param prior_rhs A named vector with the degrees of freedom and scale of the Student's t slab for regularizing large eigenvector coefficients and the prior scale for the global shrinkage parameter; e.g. \code{prior_rhs = c(slab_df = 15, slab_scale = 5, scale_global = .5)}. If \code{p0} is provided and prior_rhs is missing, \code{p0} will be used automatically to calculate \code{scale_global}.
#' @param prior_nu Set the parameters for the Gamma prior distribution on the degrees of freedom in the likelihood function when using \code{family = student_t}. Defaults to \code{prior_nu = c(alpha = 2, beta = .1)}.
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts (\code{alpha_re}). The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 20, location = 0, scale = 20)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{TRUE} for computational efficiency. This alters the interpretation of the intercept term! See \code{Details}) below.
#' @param scalex Should the covariates be scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 5000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model. Parameters from the RHS prios include \code{tau} (the global shrinkage parameter) and \code{lambda} (the local shrinkage parameter).
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} excep that \code{adapt_delta} is raised to \code{.99} and \code{max_treedepth = 15}.
#' @param zero.policy For \link[spdep]{poly2nb} which is used internally to create a spatial weights matrix. Default \code{zero.policy = TRUE}
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details If the \code{centerx = TRUE} (the default), then the intercept is the expected value of the outcome variable when 
#'   all of the covariates are at their mean value. This often has interpretive value in itself
#'   though it is the  default here for computational reasons. 
#'    
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
#'  mean log pointwise predictive density (\code{lpd}), residual spatial autocorrelation (Moran coefficient of the residuals), 
#'   root mean square error (RMSE), and median absolute deviation of residuals. Residuals are taken at the median value for each observation.}
#' \item{data}{a data frame containing the model data}
#' \item{EV}{A matrix of eigenvectors created with \code{w} and \code{geostan::make_EV}}
#' \item{C}{The spatial weights matrix used to construct EV}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re},  the random effects (varying intercepts) formula if provided, and 
#'  \code{data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("esf") and method ("ESF")}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source 
#'
#' Chun, Y., Griffith, D. A., Lee, M., and Sinha, P. (2016). "Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters." Journal of Geographical Systems, 18(1), 67-85. DOI: 10.1007/s10109-015-0225-3
#'
#' Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196(3-4), 483-493.
#' 
#' Donegan, C. (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. OSF Preprints. January 18. https://osf.io/fah3z
#'
#' Griffith, Daniel A., and Peres-Neto, Pedro R. (2006). "Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses." Ecology 87(10), 2603-2613.
#' 
#' Griffith, D., and Chun, Y. (2014). Spatial autocorrelation and spatial filtering, Handbook of Regional Science. Fischer, MM and Nijkamp, P. eds.
#'
#' Piironen, J and Vehtari, A. (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In Electronic Journal of Statistics, 11(2):5018-5051. https://projecteuclid.org/euclid.ejs/1513306866
#' 
#' @examples 
#' \dontrun{
#' library(rstan)
#' library(bayesplot)
#' library(ggplot2)
#' library(sf)
#' options(mc.cores = parallel::detectCores())
#' 
#' ## model 2016 Presidential election results in Ohio.
#' ## models here have 1 chain and small iter, only for speed of compilation. 
#' data(turmoil)
#' ohio <- turmoil[which(turmoil$state_po == "OH"),]
#' fit <- stan_esf(gop_growth ~ historic_gop + log(population) + college_educated,
#'                 slx = ~ historic_gop + college_educated,
#'                 shape = ohio,
#'                 chains = 1, iter = 1e3, family = student_t())
#' 
#' ## trace plots
#' plot(fit, plotfun = 'trace')
#' 
#' ## Rhat diagnostics (should all be very near 1; greater than 1.1 is bad).
#'   # this will include all parameters and generated quantities (residuals, log likelihood, etc.)
#' stan_rhat(fit$stanfit)
#' 
#' ## to see the eigenvector coefficients only
#' stan_rhat(fit$stanfit, par = "beta_ev")
#' 
#' ## check effective sample size of the posterior draws
#' stan_ess(fit$stanfit)
#' 
#'## summary of estimates
#' fit$summary
#' print(fit)
#' plot(fit, pars = c("beta", "intercept", "sigma"))
#' 
#' ## map the spatial filter
#' ohio$esf <- spatial(fit)$mean
#' plot(ohio[,'esf'])
#' 
#' ## Moran plot of residuals (looking for residual spatial autocorrelation)
#' w <- shape2mat(ohio, "W")
#' res_mc <- resid(fit)$mean
#' moran_plot(res_mc, w)
#' 
#' ## set priors for the main coefficients:
#'  ## coefficient priors must be a matrix with 3 columns (degrees of freedom, location, scale)
#' df <- rep(15, times = 3)
#' loc <- rep(0, times=3)
#' scale <- rep(1,times=3)
#' priors <- cbind(df, loc, scale)
#' ## with scalex=TRUE all the covariates will be centered and scaled to have stand. dev = 1
#' fit <- stan_esf(gop_growth ~ historic_gop + log(population) + college_educated, 
#'                prior = priors,
#'                shape = ohio, iter = 3e3, chains = 4, scalex = TRUE)
#'
#' ## Poisson models. Model Jim Crow era prison sentencing risk in Florida.
#' data(sentencing)
#' fit <- stan_esf(sents ~ offset(expected_sents), re = ~ name, family = poisson(),
#'                 shape = sentencing, chains = 1, iter = 1e3)
#' plot(fit, plotfun = 'trace')
#' stan_rhat(fit$stanfit, "beta_ev") +
#'    theme_bw()
#' 
#' ## map relative risk ratios
#' E <- sentencing$expected_sents
#' sentencing$ssr <- fitted(fit)$mean / E
#' ggplot(st_as_sf(sentencing),
#'        aes(fill = ssr)) +
#'        geom_sf() +
#'        scale_fill_gradient2(midpoint = 1) +
#'        theme_bw() +
#'        ggtitle("Relative risk of sentencing, 1905-1910")
#' }
#'
stan_esf <- function(formula, slx, re, shape, EV, data, C, nsa = FALSE, threshold = 0.2, family = gaussian(), p0,
                     prior = NULL, prior_intercept = NULL, prior_sigma = NULL, prior_rhs = NULL, prior_nu = NULL,
                     prior_tau = NULL,
                centerx = TRUE, scalex = FALSE, chains = 4, iter = 5e3, refresh = 500, pars = NULL,
                control = list(adapt_delta = .99, max_treedepth = 15), zero.policy = TRUE, ...) {
  if (class(family) != "family" | !family$family %in% c("gaussian", "student_t", "poisson")) stop ("Must provide a valid family object: gaussian(), student_t(), or poisson().")
  if (missing(formula)) stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(shape) & ( missing(EV) | missing(data)) ) stop ("Must provide spatially referenced data (shape) or a matrix of eigenvectors (EV) and data.")
  if(!missing(shape)) {
    shape_class <- class(shape)
    if (!any(c("sf", "SpatialPolygonsDataFrame") %in% shape_class)) stop ("shape must be of class sf or SpatialPolygonsDataFrame.")
    if ("SpatialPolygonsDataFrame" %in% shape_class) tmpdf <- shape@data
    if ("sf" %in% shape_class) tmpdf <- as.data.frame(shape)
    C <- shape2mat(shape, style = "B")  
    EV <- make_EV(C, nsa = nsa, threshold = threshold)
    dev <- ncol(EV)
    n <- nrow(EV)
  } else {
      tmpdf <- as.data.frame(data)
      dev <- ncol(EV)
      n <- nrow(tmpdf)
    }
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
        if (!exists("C")) stop("To create slx term, you must provide either C (a connectivity matrix) or shape (a Simple Features sf object or a SpatialPolygonsDataFrame).")
           Wx <- SLX(f = slx, DF = tmpdf, SWM = C, cx = centerx, sx = scalex)
           x <- cbind(Wx, x)
    } 
    dx <- ncol(x)
  }
  ModData <- make_data(formula, tmpdf, x)
  frame <- model.frame(formula, tmpdf)
  y <- model.response(frame)
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
    has_re <- 1
    id <- tmpdf[,paste(re[2])]
    n_ids <- length(unique(id))
    id_index <- to_index(id)
    re_list <- list(formula = re, data = id_index)
  } 
  if (family$family %in% c("poisson")) rhs_scale_global = 1
  if (family$family %in% c("gaussian", "student_t")) {
      # if no prior is provided at all for the rhs:
      if(missing(p0) & is.null(prior_rhs)) {
          if (exists("C")) {
              p0 <- exp_pars(formula = formula, data = tmpdf, C = C)
              if(p0 >= ncol(EV)) p0 <- ncol(EV) - .1
              rhs_scale_global <- min(1, p0/(dev-p0) / sqrt(n))
          } else{ # possible to fail without the message due to evaluation of C
              stop("p0 and prior_rhs are both missing; without providing a spatial weights matrix C, the prior for this model cannot be set.")
        }
     }    
  }
  if (!is.null(prior_rhs)) rhs_scale_global <- as.numeric(prior_rhs["scale_global"])
  is_student <- family$family == "student_t"
  priors <- list(intercept = prior_intercept, beta = prior, sigma = prior_sigma, nu = prior_nu, rhs = prior_rhs, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = priors, y = y, x = x, xcentered = centerx,
                         rhs_scale_global = rhs_scale_global, link = family$link, EV = EV)
  standata <- list(
    y = y,
    x = x,
    EV = EV,
    n = n,
    dev = dev,
    dx = dx,
    dim_beta_prior = max(1, dx),
    offset = offset,
    has_re = has_re,
    n_ids = n_ids,
    id = id_index$idx,
    slab_scale = priors$rhs["slab_scale"],
    slab_df = priors$rhs["slab_df"],
    scale_global = priors$rhs["scale_global"],
    alpha_prior = priors$intercept,
    beta_prior = t(priors$beta), 
    sigma_prior = priors$sigma,
    alpha_tau_prior = priors$alpha_tau,
    t_nu_prior = priors$nu,
    is_student = is_student
  )
  pars <- c(pars, 'intercept', 'esf', 'residual', 'beta_ev', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
  if (is_student) pars <- c(pars, "nu")
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  priors <- priors[which(names(priors) %in% c(pars, "rhs"))]
  if (family$family %in% c("gaussian", "student_t")) {
    if (intercept_only) {
      samples <- rstan::sampling(stanmodels$esf_io, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
     } else {
       samples <- rstan::sampling(stanmodels$esf_continuous, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
     }
    }
  if (family$family == "poisson") {
      samples <- rstan::sampling(stanmodels$esf_count, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, init_r = 1, ...)
  }
  out <- clean_results(samples, pars, is_student, has_re, C, x)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$C <- C
  out$EV <- EV
  out$re <- re_list
  out$priors <- priors
  out$spatial <- data.frame(par = "esf", method = "RHS-ESF")
  class(out) <- append("geostan_fit", class(out))
  return (out)
}



