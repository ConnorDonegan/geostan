#' Simultaneous autoregressive (SAR) models
#'
#' @description Fit data to a simultaneous spatial autoregressive (SAR) model, or use the SAR model as the prior model for a parameter vector in a hierarchical model.
#'
#' @param formula A model formula, following the R \code{\link[stats]{formula}} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' 
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms. 
#' 
#' @param re To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:
#' ```
#' alpha_re ~ N(0, alpha_tau)
#' alpha_tau ~ Student_t(d.f., location, scale).
#' ```
#' With the SAR model, any \code{alpha_re} term should be at a *different* level or scale than the observations; that is, at a different scale than the autocorrelation structure of the SAR model itself.
#'
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'
#' @param C Spatial weights matrix (conventionally referred to as \eqn{W} in the SAR model). Typically, this will be created using `geostan::shape2mat(shape, style = "W")`. This will be passed internally to \code{\link[geostan]{prep_sar_data}}, and will also be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms. See \code{\link[geostan]{shape2mat}}.
#'
#' @param sar_parts List of data constructed by \code{\link[geostan]{prep_sar_data}}. If not provided, then `C` will automatically be passed to \code{\link[geostan]{prep_sar_data}} to create `sar_parts`.
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{auto_gaussian()}, \code{binomial()} (with logit link function) and \code{poisson()} (with log link function); if `family = gaussian()` is provided, it will automatically be converted to `auto_gaussian()`.
#'
#' @param type Type of SAR model (character string): spatial error model ('SEM'), spatial Durbin error model ('SDEM'), spatial Durbin lag model ('SDLM'), or spatial lag model ('SLM'). see Details below.
#'
#' @param prior A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):
#' \describe{
#' 
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model `formula`. Note that if you also use `slx` terms (spatially lagged covariates), and you use custom priors for `beta`, then you have to provide priors for the slx terms. Since slx terms are *prepended* to the design matrix, the prior for the slx term will be listed first.
#' }
#'
#' \item{sar_scale}{Scale parameter for the SAR model, \code{sar_scale}. The scale is assigned a Student's t prior model (constrained to be positive).}
#'
#' \item{sar_rho}{The spatial autocorrelation parameter in the SAR model, `rho`, is assigned a uniform prior distribution. By default, the prior will be uniform over all permissible values as determined by the eigenvalues of the spatial weights matrix. The range of permissible values for `rho` is printed to the console by \code{\link[geostan]{prep_sar_data}}.}
#'
#' \item{tau}{The scale parameter for any varying intercepts (a.k.a exchangeable random effects, or partial pooling) terms. This scale parameter, `tau`, is assigned a Student's t prior (constrained to be positive).}
#' 
#' }
#'
#' @param ME To model observational uncertainty in any or all of the covariates (i.e. measurement or sampling error), provide a list of data constructed by the \code{\link[geostan]{prep_me_data}} function. 
#'
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. This increases sampling speed. If the ME argument is used, the modeled covariate (i.e., the latent variable), rather than the raw observations, will be centered. 
#'
#' @param censor_point Integer value indicating the maximum censored value; this argument is for modeling censored (suppressed) outcome data, typically disease case counts or deaths which are left-censored to protect confidentiality when case counts are very low.
#'
#' @param zmp Use zero-mean parameterization for the SAR model? Only relevant for Poisson and binomial outcome models (i.e., hierarchical models). See details below; this can sometimes improve MCMC sampling when the data is sparse, but does not alter the model specification.
#' 
#' @param prior_only Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.
#' @param chains Number of MCMC chains to use. 
#' @param iter Number of MCMC samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.
#' @param pars Specify any additional parameters you'd like stored from the Stan model.
#' @param keep_all  If `keep_all = TRUE` then samples for all parameters in the Stan model will be kept; this is necessary if you want to do model comparison with Bayes factors using the `bridgesampling` package.
#' @param slim If `slim = TRUE`, then the Stan model will not save the most memory-intensive parameters (including n-length vectors of fitted values, other 'random effects', and ME-modeled covariate values). This will disable some convenience functions that are otherwise available for fitted \code{geostan} models, such as the extraction of residuals, fitted values, and spatial trends, spatial diagnostics, and ME diagnostics. The "slim" option is designed for data-intensive routines, such as regression with raster data, Monte Carlo studies, and measurement error models. 
#' @param drop Provide a vector of character strings to specify the names of any parameters that you do not want MCMC samples for. Dropping parameters in this way can improve sampling speed and reduce memory usage. The following parameter vectors can potentially be dropped from SAR models:
#' \describe{
#' \item{fitted}{The N-length vector of fitted values}
#' \item{alpha_re}{Vector of 'random effects'/varying intercepts.}
#' \item{log_lambda_mu}{Linear predictor inside the SAR model (for Poisson and binomial models)}
#' \item{x_true}{N-length vector of 'latent'/modeled covariate values created for measurement error (ME) models.}
#' }
#' Using `drop = c('fitted', 'alpha_re', 'x_true', 'log_lambda_mu')` is equivalent to `slim = TRUE`. Note that if `slim = TRUE`, then `drop` will be ignored---so only use one or the other.
#' @param control A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details. 
#' 
#' @param ... Other arguments passed to \code{\link[rstan]{sampling}}.
#'
#' @param quiet Controls (most) automatic printing to the console. By default, any prior distributions that have not been assigned by the user are printed to the console; if `quiet = TRUE`, these will not be printed. Using `quiet = TRUE` will also force `refresh = 0`.
#' 
#' @details
#'
#' Discussions of SAR models may be found in Cliff and Ord (1981), Cressie (2015, Ch. 6), LeSage and Pace (2009), and LeSage (2014). The Stan implementation draws from Donegan (2021).
#'
#' There are two SAR specification options which are commonly known as the spatial error ('SEM') and the spatial lag ('SLM') models. When the spatial-lag of the covariates are included, then the model is referred to as a spatial Durbin model; depending on the model type, it becomes a spatial Durbin error model ('SDEM') or a spatial Durbin lag model ('SDLM').
#'
#' The spatial error specification ('SEM') is
#' \deqn{y = \mu + ( I - \rho C)^{-1} \epsilon}
#' \deqn{\epsilon \sim Gauss(0, \sigma^2)}
#' where \eqn{C} is the spatial connectivity matrix, \eqn{I} is the n-by-n identity matrix, and \eqn{\rho} is a spatial autocorrelation parameter. In words, the errors of the regression equation are spatially autocorrelated.
#'
#' Re-arranging terms, the model can also be written as follows:
#' \deqn{y = \mu + \rho C (y - \mu)  + \epsilon}
#' which perhaps shows more intuitively the implicit spatial trend component, \eqn{\rho C (y - \mu)}. 
#'
#' The second SAR specification type is the 'spatial lag of y' ('SLM'). This model describes a diffusion or contagion process:
#' \deqn{y = \rho C y + \mu + \epsilon}
#' \deqn{\epsilon \sim Gauss(0, \sigma^2)}
#' This is attractive for modeling actual contagion processes. Here the 'spatial trend' part is simply \eqn{\rho C y}.
#'
#' Both model types have a covariance matrix of:
#' 
#' \deqn{\Sigma = \sigma^2 (I - \rho C)^{-1}(I - \rho C')^{-1}.}
#'
#' But the expected values of the models differ. The expected value for the SEM is the usual \eqn{\mu} (the intercept plus \code{X*beta}); the expected value of the SLM is \eqn{(I - \rho C)^{-1} \mu}. 
#' 
#' The mathematics and typical interpretation of the SLM/SDLM is unusual and the conventional interpretation of regression coefficients does not apply! Use the \link[geostan]{impacts} method to interpret results from the SLM and SDLM models (that is, granted that this model form is at least plausible for the application).
#'
#' Use of the 'SDEM' and 'SDLM' options are for convenience: you can also obtain the Durbin models using the \code{slx} (spatial lag of X) argument. The \code{slx} argument allows control over which covariates will be added in spatial-lag form; the Durbin options include the spatial lag of all covariates.
#'
#' Most often, the SAR model is applied directly to observations (referred to below as the auto-normal or auto-Gaussian model). The SAR model can also be applied to a vector of parameters inside a hierarchical model. The latter enables spatial or network autocorrelation to be modeled when the observations are discrete counts (e.g., hierarchical models for disease incidence rates). Currently these hierarchical models are only supported for the spatial error models (SEM/SDEM).
#' 
#' ###  Auto-normal: spatial error
#'
#' When \code{family = auto_gaussian()} and `type = 'SEM'` (the default), the SAR model is specified as follows:
#' 
#' \deqn{y \sim Gauss(\mu, \Sigma)}
#' \deqn{\Sigma = \sigma^2 (I - \rho C)^{-1}(I - \rho C')^{-1}}
# 
#' where \eqn{\mu} is the mean vector (with intercept, covariates, etc.), \eqn{C} is a spatial weights or connectivity matrix (usually row-standardized), and \eqn{\sigma} is a scale parameter.
#'
#' The SAR model contains an implicit spatial trend (i.e., spatial autocorrelation) component \eqn{\phi} which is calculated as follows:
#' \deqn{
#' \phi = \rho C (y - \mu)
#' }
#' 
#' This term can be extracted from a fitted auto-Gaussian model using the \link[geostan]{spatial} method.
#'
#' When applied to a fitted auto-Gaussian model, the \link[geostan]{residuals.geostan_fit} method returns 'de-trended' residuals \eqn{R} by default. That is,
#' \deqn{
#' R = y - \mu - \rho C (y - \mu).
#' }
#' To obtain "raw" residuals (\eqn{y - \mu}), use `residuals(fit, detrend = FALSE)`. Similarly, the fitted values obtained from the \link[geostan]{fitted.geostan_fit} will include the spatial trend term by default.
#'
#' ### Aut-normal: spatial lag
#'
#'  For options `type = 'SLM'` and `type = 'SDLM'`, the \link[geostan]{spatial} method returns the vector
#' \deqn{ \phi = \rho C y, }
#' the spatial lag of \eqn{y}.
#'
#' The \link[geostan]{residuals.geostan_fit} method returns 'de-trended' residuals \eqn{R} by default:
#' \deqn{R = y - \rho C y - \mu,}
#' where \eqn{\mu} contains the intercept and any covariates (and possibly other terms). Similarly, the fitted values obtained from the \link[geostan]{fitted.geostan_fit} will include the spatial trend \eqn{\rho C y} by default. 
#'
#' To read/interpret results from the SLM or SDLM, use the \link[geostan]{impacts} method.
#' 
#' ### Poisson
#'
#' For \code{family = poisson()}, the model is specified as:
#'
#' \deqn{y \sim Poisson(e^{O + \lambda})}
#' \deqn{\lambda \sim Gauss(\mu, \Sigma)}
#' \deqn{\Sigma = \sigma^2 (I - \rho C)^{-1}(I - \rho C')^{-1}.}
#' 
#' `O` is a constant/offset term. If the raw outcome consists of a rate \eqn{\frac{y}{p}} with observed counts \eqn{y} and denominator \eqn{p} (often this will be the size of the population at risk), then the offset term \eqn{O=log(p)} is the log of the denominator.
#' 
#' This same model can be written (equivalently) as:
#' 
#' \deqn{y \sim Poisson(e^{O + \mu + \phi})}
#' \deqn{ \phi \sim Gauss(0, \Sigma) }
#'
#' This second version is referred to here as the zero-mean parameterization (ZMP), since the SAR model is forced to have mean of zero. Although the non-ZMP is typically better for MCMC sampling, use of the ZMP can greatly improve MCMC sampling *when the data is sparse*. Use `zmp = TRUE` in `stan_sar` to apply this specification. (See the geostan vignette on 'custom spatial models' for full details on implementation of the ZMP.)
#' 
#' For Poisson models, the \link[geostan]{spatial} method returns the (zero-mean) parameter vector \eqn{\phi}. When `zmp = FALSE` (the default), \eqn{\phi} is obtained by subtraction: \eqn{\phi = \lambda - \mu}.
#'
#' In the Poisson SAR model, \eqn{\phi} contains a latent (smooth) spatial trend as well as additional variation around it. If you would like to extract the latent/implicit spatial trend from \eqn{\phi}, you can do so by calculating:
#' \deqn{
#'  \rho C \phi.
#' }
#' 
#' ### Binomial
#' 
#' For `family = binomial()`, the model is specified as:
#' 
#' \deqn{y \sim Binomial(N, \lambda) }
#' \deqn{logit(\lambda) \sim Gauss(\mu, \Sigma) }
#' \deqn{\Sigma = \sigma^2 (I - \rho C)^{-1}(I - \rho C')^{-1}.}
#' 
#' where outcome data \eqn{y} are counts, \eqn{N} is the number of trials, and \eqn{\lambda} is the rate of 'success'. Note that the model formula should be structured as: `cbind(sucesses, failures) ~ 1` (for an intercept-only model), such that `trials = successes + failures`.
#' 
#' For fitted Binomial models, the \code{\link[geostan]{spatial}} method will return the parameter vector \code{phi}, equivalent to:
#' 
#' \deqn{\phi = logit(\lambda) - \mu.}
#'
#' The zero-mean parameterization (ZMP) of the SAR model can also be applied here (see the Poisson model for details); ZMP provides an equivalent model specification that can improve MCMC sampling when data is sparse.
#' 
#' As is also the case for the Poisson model, \eqn{\phi} contains a latent spatial trend as well as additional variation around it. If you would like to extract the latent/implicit spatial trend from \eqn{\phi}, you can do so by calculating:
#' \deqn{
#' \rho C \phi.
#' }
#' 
#' ## Additional functionality
#'
#' The SAR models can also incorporate spatially-lagged covariates, measurement/sampling error in covariates (particularly when using small area survey estimates as covariates), missing outcome data, and censored outcomes (such as arise when a disease surveillance system suppresses data for privacy reasons). For details on these options, please see the Details section in the documentation for \link[geostan]{stan_glm}.
#' 
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame.}
#' \item{diagnostic}{Residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including CAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the varying intercepts (\code{re}) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' 
#' \item{x_center}{If covariates are centered internally (`centerx = TRUE`), then `x_center` is a numeric vector of the values on which covariates were centered.}
#' 
#' \item{spatial}{A data frame with the name of the spatial component parameter (either "phi" or, for auto Gaussian models, "trend") and method ("SAR")}
#' \item{ME}{A list indicating if the object contains an ME model; if so, the user-provided ME list is also stored here.}
#' \item{C}{Spatial weights matrix (in sparse matrix format).}
#' \item{sar_type}{Type of SAR model: 'SEM', 'SDEM', 'SDLM', or 'SLM'.}
#' }
#'
#' @author Connor Donegan, \email{connor.donegan@gmail.com}
#' 
#' @source
#'
#' Cliff, A D and Ord, J K (1981). *Spatial Processes: Models and Applications*. Pion.
#' 
#' Cressie, Noel (2015 (1993)). *Statistics for Spatial Data*. Wiley Classics, Revised Edition.
#'
#' Cressie, Noel and Wikle, Christopher (2011). *Statistics for Spatio-Temporal Data*. Wiley.
#'
#' Donegan, Connor (2021). Building spatial conditional autoregressive (CAR) models in the Stan programming language. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#'
#' LeSage, James (2014). What Regional Scientists Need to Know about Spatial Econometrics. *The Review of Regional Science* 44: 13-32 (2014 Southern Regional Science Association Fellows Address).
#' 
#' LeSage, James, & Pace, Robert Kelley (2009). *Introduction to Spatial Econometrics*. Chapman and Hall/CRC.
#'
#' @examples
#'
#' ##
#' ## simulate SAR data on a regular grid
#' ##
#' 
#' sars <- prep_sar_data2(row = 10, col = 10, quiet = TRUE)
#' w <- sars$W
#' 
#' # draw x
#' x <- sim_sar(w = w, rho = 0.5)
#'
#' # draw y = mu + rho*W*(y - mu) + epsilon
#' # beta = 0.5, rho = 0.5
#' y <- sim_sar(w = w, rho = .5, mu = 0.5 * x)
#' dat <- data.frame(y = y, x = x)
#'
#' ##
#' ## fit SEM
#' ##
#' 
#' fit_sem <- stan_sar(y ~ x, data = dat, sar = sars,
#'                     chains = 1, iter = 800)
#' print(fit_sem)
#'
#' ##
#' ## data for SDEM
#' ##
#' 
#' # mu = x*beta + wx*gamma; beta=1, gamma=-0.25
#' x <- sim_sar(w = w, rho = 0.5)
#' mu <- 1 * x - 0.25 * (w %*% x)[,1]
#' y <- sim_sar(w = w, rho = .5, mu = mu)
#' # or for SDLM:
#' # y <- sim_sar(w = w, rho = 0.5, mu = mu, type = "SLM")
#' dat <- data.frame(y=y, x=x)
#'
#' #
#' ## fit models
#' ##
#'
#' # DSEM
#' # y = mu + rho*W*(y - mu) + epsilon
#' # mu = beta*x + gamma*Wx
#' fit_sdem <- stan_sar(y ~ x, data = dat,
#'                     sar_parts = sars, type = "SDEM",
#'                     iter = 800, chains = 1,
#'                     quiet = TRUE)
#'
#' # SDLM
#' # y = rho*Wy + beta*x + gamma*Wx + epsilon
#' fit_sdlm <- stan_sar(y ~ x, data = dat,
#'                     sar_parts = sars,
#'                     type = "SDLM",
#'                     iter = 800,
#'                     chains = 1,
#'                     quiet = TRUE)
#'
#' # compare by DIC
#' dic(fit_sdem)
#' dic(fit_sdlm)
#' 
#' \donttest{
#' ##
#' ## Modeling mortality rates
#' ##
#'
#' # simple spatial regression
#' data(georgia)
#' W <- shape2mat(georgia, style = "W")
#' 
#' fit <- stan_sar(log(rate.male) ~ 1,
#'                 C = W,
#'                 data = georgia,
#'                 iter = 900
#'                 )
#' 
#' # view fitted vs. observed, etc.
#' sp_diag(fit, georgia)
#'
#' # A more appropriate model for count data:
#' # hierarchical spatial poisson model
#' fit2 <- stan_sar(deaths.male ~ offset(log(pop.at.risk.male)),
#'                 C = W,
#'                 data = georgia,
#'                 family = poisson(),
#'                 chains = 1, # for ex. speed only
#'                 iter = 900,
#'                 quiet = TRUE
#'                  )
#' 
#' # view fitted vs. observed, etc.
#' sp_diag(fit2, georgia)
#' }
#' @export
#' @md
#' @importFrom rstan extract_sparse_parts
#' 
stan_sar <- function(formula,
                     slx,
                     re,
                     data,
                     C,
                     sar_parts = prep_sar_data(C, quiet = TRUE),
                     family = auto_gaussian(),
                     type = c("SEM", "SDEM", "SDLM", "SLM"), 
                     prior = NULL,                      
                     ME = NULL,                     
                     centerx = FALSE,
                     prior_only = FALSE,
                     censor_point,
                     zmp,
                     chains = 4,
                     iter = 2e3,
                     refresh = 500,
                     keep_all = FALSE,
                     pars = NULL,
                     slim = FALSE,
                     drop = NULL,
                     control = NULL,
                     quiet = FALSE,                     
                     ...
                     ) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("gaussian", "auto_gaussian", "poisson", "binomial"))
    if (family$family == "gaussian" | family$family == "auto_gaussian") family <- auto_gaussian(type = "SAR")
    stopifnot(!missing(data))
    check_sar_parts(sar_parts)
    stopifnot(sar_parts$n == nrow(data))
    # silence?
    if (quiet) refresh <- 0
    #### SAR type ----
    type <- match.arg(type)
    Durbin <- grepl('SDEM|SDLM', type) 
    if(grepl("SLM|SDLM", type) & family$family != "auto_gaussian") stop("SLM/SDLM are only available as auto-normal models (not Poisson or binomial models).")    
    #### SAR type [stop]
    # C
    if (!missing(C)) {
        stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
        stopifnot(all(dim(C) == nrow(data)))
    } else {
        C <- sar_parts$W
    }

    # zero-mean constraint parameterization
    sar_parts$ZMP <- ifelse(missing(zmp), 0, zmp)
    if (family$family == 'auto_gaussian') sar_parts$ZMP <- 0
    # //!!//

    
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
    if (family$family %in% c("gaussian", "auto_gaussian")) {
        y <- y_tmp
        y_int <- trials <- rep(0, n)
        if ( length(y_index_list$y_mis_idx) > 0) stop("Missing values in y; missing values are not allow for CAR/SAR models with continuous (non-count) outcomes. You may want to try stan_esf.")
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
        if (!missing(slx) | Durbin == TRUE) {
            stop("You tried to provide a spatial lag of X (slx) or 'Durbin' term for an intercept only model. Did you intend to include a covariate?")
        }
        x_full <- xraw <- model.matrix(~ 0, data = mod_frame) 
        slx <- " "
        W.list <- list(w = 1, v = 1, u = 1)
        dwx <- 0
        wx_idx <- a.zero()
  } else {
      xraw <- model.matrix(formula, data = mod_frame)
      xraw <- remove_intercept(xraw)
      if (!Durbin & missing(slx)) {
          slx <- " "
          W.list <- list(w = 1, v = 1, u = 1)
          wx_idx = a.zero()
          dwx <- 0
          x_full <- xraw          
      } else {
          ## SAR: Durbin options ----
          if (Durbin == FALSE) stopifnot(inherits(slx, 'formula'))
          if (Durbin == TRUE) {
              if (!missing(slx)) {
                  stop("You provided an slx formula and also specified a Durbin model (SDEM or SDLM). To specify spatial lag of X (slx) terms manually via the 'slx' argument, the 'type' argument must be either SEM or SLM.")
              }
              slx = formula[-2]
          }

          W <- C
          if (!inherits(W, "sparseMatrix")) W <- as(W, "CsparseMatrix")
          
          xrs <- Matrix::rowSums(W)
          if (!all(xrs == 1)) W <- row_standardize(W, warn = !quiet, msg = "Row standardizing matrix C for spatial lag of X calculations.")
          
          # nb: in stan_sar, W is taken from sar_parts, not calculated here.
          Wx <- SLX(f = slx, DF = mod_frame, x = xraw, W = W, Durbin = Durbin)
          dwx <- ncol(Wx)
          wx_idx <- as.array( which(paste0("w.", colnames(xraw)) %in% colnames(Wx)), dim = dwx )
          x_full <- cbind(Wx, xraw)
          ## SAR: Durbin options [stop] ----
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
        #n = n, # getting n from sar_parts, below
        input_offset = offset,
        has_re = has_re,
        n_ids = n_ids,
        id = id_index$idx,
        center_x = centerx,    
        ## slx data -------------
        ## stan_sar: getting W list from sar_parts
       # W_w = as.array(W.list$w),
       # W_v = as.array(W.list$v),
       # W_u = as.array(W.list$u),
       # nW_w = length(W.list$w),        
        dwx = dwx,
        wx_idx = wx_idx
    )
    ## ADD MISSING/OBSERVED INDICES -------------  
    standata <- c(y_index_list, standata)    
    ## PRIORS & SAR DATA -------------  
    is_student <- FALSE ##family$family == "student_t"
    priors_made <- make_priors(user_priors = prior,
                               y = y[y_index_list$y_obs_idx],
                               trials = trials[y_index_list$y_obs_idx],
                               x = x_full,
                               link = family$link,
                               offset = offset[y_index_list$y_obs_idx])    
    if (is.null(prior$sar_scale)) {
        priors_made$sar_scale <- priors_made$sigma        
    } else {
        priors_made$sar_scale <- prior$sar_scale
    }
    if (is.null(prior$sar_rho)) {
        priors_made$sar_rho <- uniform(sar_parts$rho_min, sar_parts$rho_max)
    } else {
        priors_made$sar_rho <- prior$sar_rho
    }
    standata$sar_rho_lims = c(priors_made$sar_rho$lower, priors_made$sar_rho$upper)
    standata <- c(standata, sar_parts)
    standata <- append_priors(standata, priors_made)
    
     ## Error = 1; Lag = 2. ##
    standata$sar <- grepl("SLM|SDLM", type) + 1
                           
    ## EMPTY PLACEHOLDERS
    standata <- add_missing_parts(standata)
    ##standata <- c(standata, empty_icar_data(n), empty_car_data(), empty_esf_data(n))
    
    ## ME MODEL -------------
    me.list <- make_me_data(ME, xraw)
    standata <- c(standata, me.list)
    ## PARAMETERS TO KEEP, with SAR PARAMETERS [START] -------------            
    pars <- c(pars, 'intercept', 'sar_scale', 'sar_rho', 'fitted')    
    if (family_int < 6 & sar_parts$ZMP == 0) pars <- c(pars, 'log_lambda_mu')
    if (sar_parts$ZMP == 1) pars <- c(pars, "log_lambda")    
    if (!intercept_only) pars <- c(pars, 'beta')
    if (dwx) pars <- c(pars, 'gamma')
    if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
    if (me.list$has_me) {
        pars <- c(pars, "x_true", "mu_x_true", "sigma_x_true")
        if (me.list$spatial_me) {
            pars <- c(pars, "car_rho_x_true")
        } else {
            pars <- c(pars, "nu_x_true")
        }
    }
    if (slim == TRUE) drop <- c('fitted', 'log_lambda', 'log_lambda_mu', 'alpha_re', 'x_true')
    pars <- drop_params(pars = pars, drop_list = drop)
    priors_made_slim <- priors_made[which(names(priors_made) %in% pars)]
    ## PARAMETERS TO KEEP, with SAR PARAMETERS [STOP] -------------
    if (me.list$has_me) priors_made_slim$ME_model <- ME$prior    
    ## PRINT PRIORS  -------------
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
    out <- clean_results(samples, pars, is_student, has_re, Wx, xraw, me.list$x_me_idx)
    out$data <- ModData # data.frame(as.matrix(ModData)) # kills with massive data
    out$family <- family
    out$formula <- formula
    out$slx <- slx
    out$re <- re_list
    out$priors <- priors_made_slim
    out$x_center <- get_x_center(standata, samples)
    out$ME <- list(has_me = me.list$has_me, spatial_me = me.list$spatial_me)
    if (out$ME$has_me) out$ME <- c(out$ME, ME)
    if (family_int == 6) {
        out$spatial <- data.frame(par = "trend", method = "SAR") 
    } else {
        out$spatial <- data.frame(par = "phi", method = "SAR")
    }
    out$N <- length( y_index_list$y_obs_idx )
    out$missing <- y_index_list    
    out$sar_parts <- sar_parts
    out$C <- as(C, "sparseMatrix")
    out$sar_type <- type
    out$diagnostic <- list()
    if (any(pars == 'fitted')) {
        C <- as(C, "sparseMatrix")        
        R <- resid(out, summary = FALSE)
        rmc <- mean( apply(R, 1, mc, w = C, warn = FALSE, na.rm = TRUE) )
        out$diagnostic$Residual_MC <- rmc
    }

    return (out)
}

