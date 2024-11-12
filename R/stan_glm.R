#' Generalized linear models
#'
#' @export
#'
#' @md
#'
#' 
#' @description Fit a generalized linear model.
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models are specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' 
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms. 
#' 
#' @param re To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:
#' ```
#' alpha_re ~ N(0, alpha_tau)
#' alpha_tau ~ Student_t(d.f., location, scale).
#' ```
#' 
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' 
#' @param C Spatial connectivity matrix which will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms. See \code{\link[geostan]{shape2mat}}.
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{poisson(link = "log")}, \code{binomial(link = "logit")}, \code{student_t()}, and the default \code{gaussian()}.
#'
#' @param prior A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):
#' \describe{
#' 
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model `formula`. Note that if you also use `slx` terms (spatially lagged covariates), and you use custom priors for `beta`, then you have to provide priors for the slx terms. Since slx terms are *prepended* to the design matrix, the prior for the slx term will be listed first.
#' }
#'
#' \item{sigma}{For `family = gaussian()` and `family = student_t()` models, the scale parameter, `sigma`, is assigned a (half-) Student's t prior distribution. The half-Student's t prior for `sigma` is constrained to be positive.}
#'
#' \item{nu}{`nu` is the degrees of freedom parameter in the Student's t likelihood (only used when `family = student_t()`). `nu` is assigned a gamma prior distribution. The default prior is `prior = list(nu = gamma2(alpha = 3, beta = 0.2))`.
#' }
#'
#' \item{tau}{The scale parameter for random effects, or varying intercepts, terms. This scale parameter, `tau`, is assigned a half-Student's t prior. To set this, use, e.g., `prior = list(tau = student_t(df = 20, location = 0, scale = 20))`.}
#' }
#' 
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a list of data as constructed by the \code{\link[geostan]{prep_me_data}} function. 
#' 
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. If the ME argument is used, the modeled covariate (i.e., latent variable), rather than the raw observations, will be centered. When using the ME argument, this is the recommended method for centering the covariates.
#'
#' @param censor_point Integer value indicating the maximum censored value; this argument is for modeling censored (suppressed) outcome data, typically disease case counts or deaths. For example, the US Centers for Disease Control and Prevention censors (does not report) death counts that are nine or fewer, so if you're using CDC WONDER mortality data you could provide `censor_point = 9`.
#' 
#' @param prior_only Draw samples from the prior distributions of parameters only. 
#' @param chains Number of MCMC chains to estimate. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples; set \code{refresh=0} to silence this.
#' @param pars Specify any additional parameters you'd like stored from the Stan model.
#' @param keep_all  If `keep_all = TRUE` then samples for all parameters in the Stan model will be kept; this is required if you want to do model comparison with Bayes factors and the `bridgesampling` package.
#' @param slim If `slim = TRUE`, then the Stan model will not collect the most memory-intensive parameters (including n-length vectors of fitted values, log-likelihoods, and ME-modeled covariate values). This will disable many convenience functions that are otherwise available for fitted \code{geostan} models, such as the extraction of residuals, fitted values, and spatial trends, WAIC, and spatial diagnostics, and ME diagnostics; many quantities of interest, such as fitted values and spatial trends, can still be calculated manually using given parameter estimates. The "slim" option is designed for data-intensive routines, such as regression with raster data, Monte Carlo studies, and measurement error models. For more control over which parameters are kept or dropped, use the `drop` argument instead of `slim`.
#' @param drop Provide a vector of character strings to specify the names of any parameters that you do not want MCMC samples for. Dropping parameters in this way can improve sampling speed and reduce memory usage. The following parameter vectors can potentially be dropped from GLM models:
#' \describe{
#' \item{'fitted'}{The N-length vector of fitted values}
#' \item{'alpha_re'}{Vector of 'random effects'/varying intercepts.}
#' \item{'x_true'}{N-length vector of 'latent'/modeled covariate values created for measurement error (ME) models.}
#' }
#' Using `drop = c('fitted', 'alpha_re', 'x_true')` is equivalent to `slim = TRUE`. If `slim = TRUE`, then `drop` will be ignored.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. 
#' @param ... Other arguments passed to \link[rstan]{sampling}.
#' @param quiet Controls (most) automatic printing to the console. By default, any prior distributions that have not been assigned by the user are printed to the console. If `quiet = TRUE`, these will not be printed. Using `quiet = TRUE` will also force `refresh = 0`.
#' 
#' @details
#'
#' Fit a generalized linear model using the R formula interface. Default prior distributions are designed to be weakly informative relative to the data. Much of the functionality intended for spatial models, such as the ability to add spatially lagged covariates and observational error models, are also available in \code{stan_glm}. All of \code{geostan}'s spatial models build on top of the same Stan code used in \code{stan_glm}.
#'
#' ### Poisson models and disease mapping
#' 
#' In spatial statistics, Poisson models are often used to calculate incidence rates (mortality rates, or disease incidence rates) for administrative areas like counties or census tracts. If \eqn{y} are counts of cases, and \eqn{P} are populations at risk, then the crude rates are \eqn{y/P}. The purpose is to model risk \eqn{\eta} for which crude rates are a (noisy) indicator. Our analysis should also respect the fact that the amount of information contained in the observations \eqn{y/P} increases with \eqn{P}. Hierarchical Poisson models are often used to incorporate all of this information.
#'
#' For the Poisson model, \eqn{y} is specified as the outcome and the log of the population at risk `log(P)` needs to be provided as an offset term. For such a case, disease incidence across the collection of areas could be modeled as:
#'
#' \deqn{y \sim Poisson(e^{log(P) + \eta})}
#' \deqn{ \eta = \alpha + A}
#' \deqn{ A \sim Gauss(0, \tau)}
#' \deqn{\tau \sim Student(20, 0, 2)}
#' where \eqn{\alpha} is the mean log-risk (incidence rate) and \eqn{A} is a vector of (so-called) random effects, which enable partial pooling of information across observations. Covariates can be added to the model for the log-rates, such that \eqn{\eta = \alpha + X \beta + A}.
#'
#' Note that the denominator for the rates is specified as a log-offset to provide a consistent, formula-line interface to the model. Using the log-offest (as above) is equivalent to the following:
#' \deqn{
#' y \sim Poisson(P * e^{\eta})
#' }
#' where \eqn{P} is still the population at risk and it is multiplied by \eqn{e^{\eta}}, the incidence rate (risk). 
#'
#' ### Spatially lagged covariates (SLX)
#' 
#' The `slx` argument is a convenience function for including SLX terms. For example, 
#' \deqn{
#'  y = W X \gamma + X \beta + \epsilon
#' }
#' where \eqn{W} is a row-standardized spatial weights matrix (see \code{\link[geostan]{shape2mat}}), \eqn{WX} is the mean neighboring value of \eqn{X}, and \eqn{\gamma} is a coefficient vector. This specifies a regression with spatially lagged covariates. SLX terms can specified by providing a formula to the \code{slx} argument:
#' ```
#' stan_glm(y ~ x1 + x2, slx = ~ x1 + x2, \...),
#' ```
#' which is a shortcut for
#' ```
#' stan_glm(y ~ I(W \%*\% x1) + I(W \%*\% x2) + x1 + x2, \...)
#' ```
#' SLX terms will always be *prepended* to the design matrix, as above, which is important to know when setting prior distributions for regression coefficients.
#'
#' For measurement error (ME) models, the SLX argument is the only way to include spatially lagged covariates since the SLX term needs to be re-calculated on each iteration of the MCMC algorithm.
#' 
#' ### Measurement error (ME) models
#' 
#' The ME models are designed for surveys with spatial sampling designs, such as the American Community Survey (ACS) estimates. For a tutorial, see \code{vignette("spatial-me-models", package = "geostan")}.
#'
#' Given estimates \eqn{x}, their standard errors \eqn{s}, and the target quantity of interest (i.e., the unknown true value) \eqn{z}, the ME models have one of the the following two specifications, depending on the user input. If a spatial CAR model is specified, then:
#' 
#' \deqn{x \sim Gauss(z, s^2)}
#' \deqn{z \sim Gauss(\mu_z, \Sigma_z)}
#' \deqn{\Sigma_z = (I - \rho C)^{-1} M}
#' \deqn{\mu_z \sim Gauss(0, 100)}
#' \deqn{\tau_z \sim Student(10, 0, 40), \tau > 0}
#' \deqn{\rho_z \sim uniform(l, u)}
#'
#' where \eqn{\Sigma} specifies the covariance matrix of a spatial conditional autoregressive (CAR) model with scale parameter \eqn{\tau} (on the diagonal of \eqn{M}), autocorrelation parameter \eqn{\rho}, and \eqn{l}, \eqn{u} are the lower and upper bounds that \eqn{\rho} is permitted to take (which is determined by the extreme eigenvalues of the spatial connectivity matrix \eqn{C}). \eqn{M} contains the inverse of the row sums of \eqn{C} on its diagonal multiplied by \eqn{\tau} (following the "WCAR" specification).
#' 
#' For non-spatial ME models, the following is used instead:
#'\deqn{x \sim Gauss(z, s^2)}
#' \deqn{z \sim student_t(\nu_z, \mu_z, \sigma_z)}
#' \deqn{\nu_z \sim gamma(3, 0.2)}
#' \deqn{\mu_z \sim Gauss(0, 100)}
#' \deqn{\sigma_z \sim student(10, 0, 40)}
#' 
#' For strongly skewed variables, such as census tract poverty rates, it can be advantageous to apply a logit transformation to \eqn{z} before applying the CAR or Student-t prior model. When the `logit` argument is used, the first two lines of the model specification become:
#' \deqn{x \sim Gauss(z, s^2)}
#' \deqn{logit(z) \sim Gauss(\mu_z, \Sigma_z) }
#' and similarly for the Student t model:
#' \deqn{x \sim Gauss(z, s^2)}
#' \deqn{logit(z) \sim student(\nu_z, \mu_z, \sigma_z)}
#'
#'
#' ### Missing data
#'
#' For most geostan models, missing (NA) observations are allowed in the outcome variable. However, there cannot be any missing covariate data. Models that can handle missing data are: any Poisson or binomial model (GLM, SAR, CAR, ESF, ICAR), all GLMs and ESF models. The only models that cannot handle missing outcome data are the CAR and SAR models when the outcome is a continuous variable (auto-normal/Gaussian models). 
#'
#' When observations are missing, they will simply be ignored when calculating the likelihood in the MCMC sampling process (reflecting the absence of information). The estimated model parameters (including any covariates and spatial trend) will then be used to produce estimates or fitted values for the missing observations. The `fitted` and `posterior_predict` functions will work as normal in this case, and return values for all rows in your data.
#' 
#' ### Censored counts
#'
#' Vital statistics systems and disease surveillance programs typically suppress case counts when they are smaller than a specific threshold value. In such cases, the observation of a censored count is not the same as a missing value; instead, you are informed that the value is an integer somewhere between zero and the threshold value. For Poisson models (`family = poisson())`), you can use the `censor_point` argument to encode this information into your model. 
#'
#' Internally, `geostan` will keep the index values of each censored observation, and the index value of each of the fully observed outcome values. For all observed counts, the likelihood statement will be:
#' \deqn{
#' p(y_i | data, model) = poisson(y_i | \mu_i), 
#' }
#' as usual, where \eqn{\mu_i} may include whatever spatial terms are present in the model.
#'
#' For each censored count, the likelihood statement will equal the cumulative Poisson distribution function for values zero through the censor point:
#' \deqn{
#' p(y_i | data, model) = \sum_{m=0}^{M} Poisson( m | \mu_i),
#' }
#' where \eqn{M} is the censor point and \eqn{\mu_i} again is the fitted value for the \eqn{i^{th}} observation.
#' 
#' For example, the US Centers for Disease Control and Prevention's CDC WONDER database censors all death counts between 0 and 9. To model CDC WONDER mortality data, you could provide `censor_point = 9` and then the likelihood statement for censored counts would equal the summation of the Poisson probability mass function over each integer ranging from zero through 9 (inclusive), conditional on the fitted values (i.e., all model parameters). See Donegan (2021) for additional discussion, references, and Stan code.
#'
#' 
#' @return
#'
#' An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{C}{The spatial weights matrix, if one was provided by the user.}
#' \item{re}{A list containing \code{re}, the random effects (varying intercepts) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{x_center}{If covariates are centered internally (`centerx = TRUE`), then `x_center` is a numeric vector of the values on which covariates were centered.}
#' \item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
#' \item{spatial}{NA, slot is maintained for use in \code{geostan_fit} methods.}
#' }
#' 
#' @author Connor Donegan, \email{connor.donegan@gmail.com}
#'
#' @source
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Donegan, Connor (2021). Building spatial conditional autoregressive (CAR) models in the Stan programming language. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#' 
#' @examples
#' ##
#' ## Linear regression model
#' ##
#' 
#' N = 100
#' x <- rnorm(N)
#' y <- .5 * x + rnorm(N)
#' dat <- cbind(y, x)
#'
#' # no. of MCMC samples
#' iter = 600
#'
#' # fit model
#' fit <- stan_glm(y ~ x, data = dat, iter = iter, quiet = TRUE)
#'
#' # see results with MCMC diagnostics
#' print(fit)
#' 
#' ##
#' ## Custom prior distributions
#' ##
#'
#' PL <- list(
#'       intercept = normal(0, 1),
#'       beta = normal(0, 1),
#'       sigma = student_t(10, 0, 2)
#' )
#'
#' fit2 <- stan_glm(y ~ x, data = dat, prior = PL, iter = iter,
#'                 quiet = TRUE)
#'
#' print(fit2)
#' 
#' ##
#' ## Poisson model for count data
#' ## with county 'random effects' 
#' ##
#'
#' data(sentencing)
#'
#' # note: 'name' is county identifier
#' head(sentencing)
#' 
#' # denominator in standardized rate Y/E
#' # (observed count Y over expected count E)
#' # (use the log-denominator as the offest term)
#' sentencing$log_e <- log(sentencing$expected_sents)
#'
#' # fit model
#' fit.pois <- stan_glm(sents ~ offset(log_e),
#'                      re = ~ name,
#'                      family = poisson(),
#'                      data = sentencing,                    
#'                     iter = iter, quiet = TRUE) 
#'
#' # Spatial autocorrelation/residual diagnostics
#' sp_diag(fit.pois, sentencing)
#'
#' # summary of results with MCMC diagnostics
#' print(fit.pois)
#' 
#' # MCMC diagnostics plot: Rhat values should all by very near 1
#' rstan::stan_rhat(fit.pois$stanfit)
#' 
#' # effective sample size for all parameters and generated quantities
#' # (including residuals, predicted values, etc.)
#' rstan::stan_ess(fit.pois$stanfit)
#' 
#' # or for a particular parameter
#' rstan::stan_ess(fit.pois$stanfit, "alpha_re")
#'
#' ##
#' ## Visualize the posterior predictive distribution
#' ##
#'
#' # plot observed values and model replicate values
#' yrep <- posterior_predict(fit.pois, S = 65)
#' y <- sentencing$sents
#' ltgray <- rgb(0.3, 0.3, 0.3, 0.5)
#' 
#' plot(density(yrep[1,]), col = ltgray,
#'      ylim = c(0, 0.014), xlim = c(0, 700),
#'      bty = 'L', xlab = NA, main = NA)
#' 
#' for (i in 2:nrow(yrep)) lines(density(yrep[i,]), col = ltgray)
#' 
#' lines(density(sentencing$sents), col = "darkred", lwd = 2)
#' 
#' legend("topright", legend = c('Y-observed', 'Y-replicate'),
#'        col = c('darkred', ltgray), lwd = c(1.5, 1.5))
#'
#' # plot replicates of Y/E
#'E <- sentencing$expected_sents
#'
#' # set plot margins
#' old_pars <- par(mar=c(2.5, 3.5, 1, 1))
#'
#' # plot yrep
#' plot(density(yrep[1,] / E), col = ltgray,
#'	    ylim = c(0, 0.9), xlim = c(0, 7),
#'     bty = 'L', xlab = NA, ylab = NA, main = NA)
#'
#' for (i in 2:nrow(yrep)) lines(density(yrep[i,] / E), col = ltgray)
#'
#' # overlay y
#' lines(density(sentencing$sents / E), col = "darkred", lwd = 2)
#'
#' # legend, y-axis label
#' legend("topright", legend = c('Y-observed', 'Y-replicate'),
#'       col = c('darkred', ltgray), lwd = c(1.5, 1.5))
#'
#' mtext(side = 2, text = "Density", line = 2.5)
#'
#' # return margins to previous settings
#' par(old_pars)
#' 
#' @importFrom rstan extract_sparse_parts
stan_glm <- function(formula,
                     slx,
                     re,
                     data,
                     C,
                     family = gaussian(),
                     prior = NULL,
                     ME = NULL,                     
                     centerx = FALSE, 
                     prior_only = FALSE,
                     censor_point,
                     chains = 4,
                     iter = 2e3, 
                     refresh = 1e3,
                     keep_all = FALSE,
                     slim = FALSE,
                     drop = NULL,
                     pars = NULL,
                     control = NULL,
                     quiet = FALSE,
                     ...) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("gaussian", "student_t", "poisson", "binomial"))
    stopifnot(!missing(data))
    # silence?
    if (quiet) refresh <- 0
    # C    
    if (!missing(C)) {
        stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
        stopifnot(all(dim(C) == nrow(data)))
    } else {
        if (!missing(ME) && ME$spatial_me == TRUE) if (inherits(ME$car_parts$C, "Matrix") | inherits(ME$car_parts$C, "matrix")) C <- ME$car_parts$C
    }
    tmpdf <- as.data.frame(data)
    n <- nrow(tmpdf)    
    family_int <- family_2_int(family)        
    if (!missing(censor_point)) if (family$family != "poisson") stop("censor_point argument is only available for Poisson models.")
    if (missing(censor_point)) censor_point <- FALSE
    mod_frame <- model.frame(formula, tmpdf, na.action = NULL)
    handle_missing_x(mod_frame)
    ## prep Y for missingness //
    y_index_list <- handle_censored_y(censor_point, mod_frame)
    mi_idx <- y_index_list$y_mis_idx  
    y_tmp <- model.response(mod_frame)
    if (family$family == "binomial") {
        trials <- y_tmp[,1] + y_tmp[,2]
        y <- y_int <- y_tmp[,1]
    }
    if (family$family == "poisson") {
        y <- y_int <- y_tmp        
        trials <- rep(0, n)
    }
    if (family$family %in% c("gaussian", "student_t")) {
        y <- y_tmp
        y_int <- trials <- rep(0, n)
    }    
    y[mi_idx] <- y_int[mi_idx] <- trials[mi_idx] <- 0
    ## //
    if (is.null(model.offset(mod_frame))) {
        offset <- rep(0, times = n)
    } else {
        offset <- model.offset(mod_frame)
    }    
    x_mat_tmp <- model.matrix(formula, mod_frame)
    intercept_only <- ifelse(all(dimnames(x_mat_tmp)[[2]] == "(Intercept)"), 1, 0)
    if (intercept_only) {
        if (!missing(slx)) {
            stop("You provided a spatial lag of X (slx) term for an intercept only model. Did you intend to include a covariate?")
        }
        x_full <- xraw <- model.matrix(~ 0, data = mod_frame) 
        slx <- " "
        W.list <- list(w = 1, v = 1, u = 1)
        dwx <- 0
        wx_idx <- a.zero()
  } else {
      xraw <- model.matrix(formula, data = mod_frame)
      xraw <- remove_intercept(xraw)
      if (missing(slx)) {
          slx <- " "
          W.list <- list(w = 1, v = 1, u = 1)
          wx_idx = a.zero()
          dwx <- 0
          x_full <- xraw          
      } else {
          stopifnot(inherits(slx, "formula"))
          W <- C
          if (!inherits(W, "sparseMatrix")) W <- as(W, "CsparseMatrix")
          ## xrs <- Matrix::rowSums(W)
          ## if (!all(xrs == 1)) W <- row_standardize(W, warn = !quiet, msg = "Row standardizing matrix C for spatial lag of X calculations.")
          # efficient transform to CRS representation for W.list (via transpose)
          Wij <- as(W, "TsparseMatrix")
          Tw <- Matrix::sparseMatrix(i = Wij@j + 1,
                                     j = Wij@i + 1,
                                     x = Wij@x,
                                     dims = dim(Wij))
          W.list <- list(w = Tw@x,
                         v = Tw@i + 1,
                         u = Tw@p + 1)          
          Wx <- SLX(f = slx, DF = mod_frame, x = xraw, W = W)
          dwx <- ncol(Wx)
          wx_idx <- as.array( which(paste0("w.", colnames(xraw)) %in% colnames(Wx)), dim = dwx )
          x_full <- cbind(Wx, xraw)
      }
  }
    ModData <- make_data(mod_frame, x_full, y_index_list$y_mis_idx)
    if(missing(re)) {
        has_re <- n_ids <- id <- 0;
        id_index <- to_index(id, n = n)
        re_list <- NULL
    } else {
        stopifnot(inherits(re, "formula"))
        has_re <- 1
        id <- tmpdf[,paste(re[2])]
        n_ids <- length(unique(id))
        id_index <- to_index(id)
        re_list <- list(formula = re, data = id_index)
  }
    standata <- list(
  ## glm data -------------      
    family = family_int,
    prior_only = prior_only,    
    y = y,
    y_int = y_int,
    trials = trials,
    n = n,
    input_offset = offset,
    has_re = has_re,
    n_ids = n_ids,
    id = id_index$idx,
    center_x = centerx,      
    ## slx data -------------    
    W_w = as.array(W.list$w),
    W_v = as.array(W.list$v),
    W_u = as.array(W.list$u),
    nW_w = length(W.list$w),
    dwx = dwx,
    wx_idx = wx_idx
    )
    ## ADD MISSING/OBSERVED INDICES -------------  
    standata <- c(y_index_list, standata)
    ## PRIORS -------------
    is_student <- family$family == "student_t"    
    priors_made <- make_priors(user_priors = prior,
                               y = y[y_index_list$y_obs_idx],
                               trials = trials[y_index_list$y_obs_idx],
                               x = x_full,
                               link = family$link,
                               offset = offset[y_index_list$y_obs_idx])    
    standata <- append_priors(standata, priors_made)
    
    ## EMPTY PLACEHOLDERS
    standata <- add_missing_parts(standata)    
    ##empty_parts <- c(empty_icar_data(n), empty_esf_data(n), empty_car_data(), empty_sar_data(n))
    ##empty_parts <- empty_parts[ which(!names(empty_parts) %in% names(standata)) ]
    ##standata <- c(standata, empty_parts)
    
    ## ME MODEL -------------  
    me.list <- make_me_data(ME, xraw)
    standata <- c(standata, me.list)  
    ## PARAMETERS TO KEEP -------------               
    pars <- c(pars, 'intercept', 'fitted')
    if (!intercept_only) pars <- c(pars, 'beta')
    if (dwx) pars <- c(pars, 'gamma')
    if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
    if (is_student) pars <- c(pars, "nu")
    if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
    if (me.list$has_me) {
        pars <- c(pars, "x_true", "mu_x_true", "sigma_x_true")
        if (me.list$spatial_me) {
            pars <- c(pars, "car_rho_x_true")
        } else {
            pars <- c(pars, "nu_x_true")
        }
    }
    if (slim == TRUE) drop <- c('fitted', 'alpha_re', 'x_true')
    pars <- drop_params(pars = pars, drop_list = drop)
    priors_made_slim <- priors_made[which(names(priors_made) %in% pars)]
    if (me.list$has_me) priors_made_slim$ME_model <- ME$prior        
    if (!quiet) print_priors(prior, priors_made_slim)
    ## MCMC INITIAL VALUES ------------- 
    inits <- "random"
    if (censor_point > 0 | (standata$has_me == 1 && any(standata$use_logit == 1))) {
        FRAME <- sys.nframe()
        inits <- init_fn_builder(FRAME_NUMBER = FRAME)
    } 
    ## CALL STAN -------------
    if (keep_all == TRUE) {
        xparsx <- NA
    } else {
        xparsx <- pars
    }    
    samples <- rstan::sampling(stanmodels$foundation, data = standata, iter = iter, chains = chains, refresh = refresh, pars = xparsx, control = control, init = inits, ...) 
    ## OUTPUT -------------    
    out <- clean_results(samples, pars, is_student, has_re, Wx, xraw, me.list$x_me_idx)        
    out$data <- data.frame(as.matrix(ModData))
    out$family <- family
    out$formula <- formula
    out$slx <- slx
    out$re <- re_list
    out$priors <- priors_made_slim
    out$x_center <- get_x_center(standata, samples)
    out$ME <- list(has_me = me.list$has_me, spatial_me = me.list$spatial_me)    
    if (out$ME$has_me) out$ME <- c(out$ME, ME)
    if (has_re) {
        out$spatial <- data.frame(par = "alpha_re", method = "Exchangeable")
    } else {
        out$spatial <- data.frame(par = "none", method = "none")
    }
    if (!missing(C)) out$C <- as(C, "sparseMatrix")    
    out$N <- length( y_index_list$y_obs_idx )
    out$missing <- y_index_list

    out$diagnostic <- list()    
    if (!missing(C) && any(pars == 'fitted')) {
        C <- as(C, "sparseMatrix")        
        R <- resid(out, summary = FALSE)
        rmc <- mean( apply(R, 1, mc, w = C, warn = FALSE, na.rm = TRUE) )
        out$diagnostic$Residual_MC <- rmc
    }    

    return(out)
    
}

