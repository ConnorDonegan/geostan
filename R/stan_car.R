#' Conditional autoregressive (CAR) models
#'
#' @export
#' 
#' @description Apply the CAR model as a prior on parameters, or fit an auto Gaussian model.
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized version of the user-provided spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term \code{alpha_re} specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re ~ N(0, alpha_tau)}, \code{alpha_tau ~ Student_t(d.f., location, scale)}. With the CAR model, any \code{alpha_re} term should be at a different level or scale than the observations (i.e., a different scale than the CAR model). In other words, do not try to replicate the BYM model using the CAR model because that would be redundant.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list. Errors are assigned a Gaussian probability distribution and the modeled (true) covariate vector is assigned an auto-Gaussian (CAR) model unless \code{spatial = FALSE}, in which case they are assinged the Student's t model. Elements of the list \code{ME} may include:
#' \describe{
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one hundred) provide a vector of zeros and ones indicating which columns are bounded. By default, if \code{bounded = TRUE}, the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of any bounded variables.}
#' \item{spatial}{Logical value indicating whether an auto-Gaussian (conditional autoregressive (CAR)) model should be used for the covariates. For \code{stan_car}, defaults to \code{spatial = TRUE}.}
#' \item{car_parts}{This argument is ignored for \code{stan_car}, because the user must always provide this information within the main function call. See the \code{car_parts} argument below. If any information is provided to \code{ME$car_parts}, it will be overwritten internally.}
#' }
#' 
#' @param car_parts A list of data for the CAR model, as returned by \link[geostan]{prep_car_data} (be sure to return the matrix \code{C}, by using the argument \code{cmat = TRUE}).
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{gaussian()}, \code{binomial(link = "logit")}, and \code{poisson(link = "log")}.
#' @param invert To calculate the log likelihood of the data \code{log_lik} and the posterior predictive distribution \code{yrep} with the auto-Gaussian model, the precision matrix needs to be inverted. This can be costly for large data sets---the inversion needs to be completed once per posterior sample. To avoid the computational cost, set \code{invert = FALSE}. Note, this is only used when \code{family = gaussian()}.
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the regression coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of varying intercepts \code{alpha_re}. The varying intercepts are given a normal prior with scale parameter \code{alpha_tau}. \code{alpha_tau} is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 15, location = 0, scale = 5)}.
#' @param prior_car_scale Hyperprior parameters for the scale of the CAR model \code{car_scale}. The scale is assigned a Student's t prior model; to set its parameter values provide a length-three vector with the degrees of freedom, location, and scale parameters. E.g., \code{prior_car_scale = c(df = 15, location = 0, scale = 2)}.
#' @param centerx Logical value indicating if the covariates should be centered prior to fitting the model.
#' @param scalex Logical value indicating if the covariates be centered and scaled (divided by their standard deviation).
#' @param prior_only Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.
#' @param chains Number of MCMC chains to use. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} except that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param silent If \code{TRUE}, suppress printed messages including prior specifications and Stan sampling progress (i.e. \code{refresh=0}). Stan's error and warning messages will still print.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first. 
#' @details
#'  
#' The CAR model is specified as follows:
#'
#'                             \code{Y ~ MVGauss(Mu, (I - rho C)^-1 * M * tau^2}
#'
#' where \code{Y} may be data or parameters and \code{Mu} is the mean vector. When \code{Y} is a vector of parameters, \code{Mu} will be set to a vector of zeroes. \code{C} is the spatial connectivity matrix and \code{M} is a known diagonal matrix with diagonal entries proportional to the conditional variances (see the \code{M_diagonal} argument). There are two parameters to estimate: \code{car_rho} (rho), which controls the degree of spatial autocorrelation, and the scale parameter, \code{car_scale} (tau). 
#' 
#' These models are discussed in Cressie and Wikle (2011, p. 184-88), Cressie (2015 [1993], Ch. 6-7), and Haining and Li (2020, p. 249-51). The Stan code for the model draws from Joseph (2016) to calculate the log-determinant of the precision matrix (also see Haining and Li, p. 250) up to a constant of proportionality.
#'
#' The auto-Gaussian model returns a spatial \code{trend} component which is calculated as follows (Cressie 2015, p. 564):
#'
#' \code{trend = rho * C * (Y - Mu)}
#'
#' where \code{Mu} includes the intercept plus any covariates \code{X * beta} and any other terms that may be added to the linear predictor. For the Poisson and Binomial models, the spatial trend is simply equal to the parameter vector \code{phi} (however, note that a spatial trend could be extracted from \code{phi} in this same manner). All models return:
#'
#' \code{fitted = Mu}.
#'
#' The residuals returned by the auto-Gaussian model are `de-trended', as follows:
#' 
#' \code{residual = Y - Mu - trend}
#'
#' The trend, fitted values, and residuals may be obtained using the \link[geostan]{spatial}, \link[geostan]{fitted.geostan_fit}, and \link[geostan]{resid.geostan_fit} methods, as usual. For Poisson and Binomial models, the fitted values and residuals are calculated in the usual manner.
#'
#' The posterior predictive distribution (ppd) is returned in the parameter vector \code{yrep}. For the auto-Gaussian model, each sample from the ppd is a draw from the multivariate Gaussian density function parameterized as:
#'
#' \code{yrep \~ MVGauss(Mu, S) },   \code{S = (I - rho C)^-1 * M * tau^2}
#'
#' and is accessible using the \link[geostan]{posterior_predict} method, as usual. The log likelihood of the data conditional on the model (\code{log_lik}, required for calculating \link[geostan]{waic}) is calculated analogously, using the above parameterization of the multivariate normal density. 
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame.}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including CAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the varying intercepts (\code{re}) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{spatial}{A data frame with the name of the spatial component parameter (either "phi" or, for auto Gaussian models, "trend") and method ("CAR")}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source
#'
#' Cressie, Noel (2015 [1993]). Statistics for Spatial Data. Wiley Classics, Revised Edition.
#' 
#' Cressie, Noel and Wikle, Christopher (2011). Statistics for Spatio-Temporal Data. Wiley.
#'
#' Haining, Robert and Li, Guangquan (2020). Modelling Spatial and Spatial-Temporal Data: A Bayesian Approach. CRC Press.
#' 
#' Joseph, Max (2016). Exact Sparse CAR Models in Stan. Stan Case Studies, Vol. 3. \url{https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html}
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(bayesplot)
#' library(sf)
#' options(mc.cores = parallel::detectCores())
#' data(sentencing)
#'
#' C <- shape2mat(sentencing, style = "B")
#' car.dl <- prep_car_data(C, style = "ACAR")
#' log_e <- log(sentencing$expected_sents)
#' fit.car <- stan_car(sents ~ offset(log_e),
#'                     family = poisson(),
#'                     data = sentencing,
#'                     car_parts = car.dl,
#'                     refresh = 0
#' )
#' 
#' # posterior predictive distribution
#' yrep <- posterior_predict(fit.car, samples = 75)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep)
#'
#' plot(fit.car, pars = c("car_scale", "car_alpha"))
#' 
#' sp.trend <- spatial(fit.car)$mean
#' st_as_sf(sentencing) %>%
#'   ggplot() +
#'   geom_sf(aes(fill = sp.trend)) +
#'   scale_fill_gradient2()
#'  }
#' 
stan_car <- function(formula,
                     slx,
                     re,
                     data,
                     ME = NULL,
                     car_parts,
                     family = gaussian(), 
                     invert = TRUE, #!#
                     prior = NULL, prior_intercept = NULL, prior_tau = NULL, prior_car_scale = NULL,
                     centerx = FALSE, scalex = FALSE,
                     prior_only = FALSE,
                     chains = 5, iter = 2e3, refresh = 500, pars = NULL,
                     control = list(adapt_delta = 0.9, max_treedepth = 13),
                     silent = FALSE,
                     ...
                     ) {
  if (!inherits(family, "family") | !family$family %in% c("binomial", "poisson", "gaussian")) stop ("Must provide a valid family object: gaussian(), binomial() or poisson().") #!#
  if (missing(formula) | !inherits(formula, "formula")) stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(data)) stop("Must provide data (a data.frame or object coercible to a data.frame).")
  if (missing(car_parts) | !inherits(car_parts, "list")) stop("Must provide a list to the argument car_parts.")
  if(!all(c("nC", "nImC", "ImC", "ImC_v", "ImC_u", "Cidx", "M_diag", "C") %in% names(car_parts))) stop("car_parts is missing at least one required part. See ?prep_car_data. Did you use cmat = TRUE?")
  C <- car_parts$C
  if (!inherits(C, "matrix")) stop("car_parts$C must be a matrix.")
  if (scalex) centerx = TRUE
  if (silent) refresh = 0
  ## GLM STUFF -------------
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")    
  n <- nrow(mod.mat)
  if (any(dim(C) != n)) stop("Dimensions of matrix C must match the number of observations. See ?shape2mat and ?prep_car_data for help creating C.")      
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
  is_student <- FALSE # family$family == "student_t" # always false
  user_priors <- list(intercept = prior_intercept,
                      beta = prior,
                      sigma = NULL,
                      nu = NULL,
                      alpha_tau = prior_tau)
  priors <- make_priors(user_priors = user_priors,
                        y = y,
                        x = x,
                        xcentered = centerx,
                        link = family$link,
                        offset = offset)
  ## CAR STUFF -------------  #!#
  if (family$family == "gaussian") family_int <- 5 #!# auto-Gaussian #!#
  if (!is.null(prior_car_scale)) {
      priors$sigma <- prior_car_scale
  } else {
      prior_car_scale <- priors$sigma
      if (!silent) {
          message("\n*Setting prior parameters for car_scale")
          message("Student's t")
          message("Degrees of freedom: ", priors$sigma[1])
          message("Location: ", priors$sigma[2])
          message("Scale: ", priors$sigma[3])
      }      
  }  
  ## MIXED STUFF -------------  
  standata <- list(
  ## glm data -------------
    y = y,
    y_int = y_int,
    trials = rep(0, length(y)),
    n = n,
    dbeta_prior = dbeta_prior,
    offset = offset,
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
  ## draw samples from prior only, ignoring data and likelihood
    prior_only = prior_only
  )
  ## CAR STUFF
  standata <- c(standata, car_parts)
  ## DATA MODEL STUFF -------------
   ## overwrites any use specified ME$car_parts
  if (!is.null(ME)) ME$car_parts <- car_parts
  me.list <- prep_me_data(ME, x.list$x)
  standata <- c(standata, me.list)
  ## STAN STUFF -------------    
  if (family$family == "binomial") {
      # standata$y will be ignored for binomial and poisson models
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'car_scale', 'car_rho', 'residual', 'fitted')
  if (invert) pars <- c(pars, 'log_lik', 'yrep')
  if (family_int < 5) pars <- c(pars, 'phi') #!#
  if (family_int == 5) pars <- c(pars, "trend")
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  priors <- priors[which(names(priors) %in% pars)]
  ## CAR SCALE PRIOR -------------    
  priors <- c(priors, list(car_scale = prior_car_scale)) #!#
  ## PRINT STUFF -------------    
  if (!silent) print_priors(user_priors, priors)
  ## CALL STAN -------------  
  samples <- rstan::sampling(stanmodels$car, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  out <- clean_results(samples, pars, is_student, has_re, car_parts$C, Wx, x.list$x, me.list$x_me_unbounded_idx, me.list$x_me_bounded_idx)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME
  if (family_int != 5) out$spatial <- data.frame(par = "phi", method = "CAR") else out$spatial <- data.frame(par = "trend", method = "CAR") #!#
  class(out) <- append("geostan_fit", class(out))
  return (out)
}

