#' Conditional autoregressive (CAR) models
#' 
#' @description Use the CAR model as a prior on parameters, or fit data to an auto-Gaussian CAR model.
#' 
#' @param formula A model formula, following the R \code{\link[stats]{formula}} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' 
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms. 
#' 
#' @param re To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:
#' ```
#'        alpha_re ~ N(0, alpha_tau)
#'        alpha_tau ~ Student_t(d.f., location, scale).
#' ```
#' With the CAR model, any \code{alpha_re} term should be at a *different* level or scale than the observations; that is, at a different scale than the autocorrelation structure of the CAR model itself.
#'
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' 
#' @param car_parts A list of data for the CAR model, as returned by \code{\link[geostan]{prep_car_data}}.
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{auto_gaussian()}, \code{binomial(link = "logit")}, and \code{poisson(link = "log")}; if `family = gaussian()` is provided, it will automatically be converted to `auto_gaussian()`.
#'
#' @param prior A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):
#' \describe{
#' 
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model `formula`. Note that if you also use `slx` terms (spatially lagged covariates), and you use custom priors for `beta`, then you have to provide priors for the slx terms. Since slx terms are *prepended* to the design matrix, the prior for the slx term will be listed first.
#' }
#'
#' \item{car_scale}{The scale of the CAR model, \code{car_scale}. The scale is assigned a Student's t prior model (constrained to be positive).}
#'
#' \item{car_rho}{The spatial autocorrelation parameter in the CAR model, `rho`, is assigned a uniform prior distribution. By default, the prior will be uniform over all permissible values as determined by the eigenvalues of the connectivity matrix, `C`. The range of permissible values for `rho` is automatically printed to the console by \code{\link[geostan]{prep_car_data}}.}
#'
#' \item{tau}{The scale parameter for any varying intercepts (a.k.a exchangeable random effects, or partial pooling) terms. This scale parameter, `tau`, is assigned a Student's t prior (constrained to be positive).}
#' 
#' }
#' 
#'@param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list of priors. See the Details section below for more information. Elements of the \code{ME} list may include:
#'
#' \describe{
#' \item{se}{A required data.frame with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' 
#' \item{bounds}{An optional numeric vector of length two providing the upper and lower bounds, respectively, of the variables. If not provided, they will be set to `c(-Inf, Inf)` (i.e., unbounded). Common usages include keeping percentages between zero and one hundred or proportions between zero and one.}
#' 
#'\item{car_parts}{For CAR models, the default is for ME models to be spatial CAR models themselves. However, to use a non-spatial Student's t ME model instead, provide `ME$car_parts = FALSE`.}
#' 
#'  \item{prior}{Optionally provide a named list of prior distributions for the measurement error model(s) (see \code{\link[geostan]{priors}}). If none are provided, default priors will be assigned and printed to the console. \describe{
#' 
#'  \item{df}{If using a non-spatial ME model, the degrees of freedom for the Student's t model is assigned a gamma prior with default parameters of `gamma(alpha = 3, beta = 0.2)`. Provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#' 
#'   \item{location}{The prior for the location parameter, mu, is Gaussian (the default being `normal(location = 0, scale = 100)`). To adjust the prior distributions, provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#' 
#'  \item{scale}{The prior for the `scale` parameters is Student's t, and the default is `student_t(df = 10, location = 0, scale = 40)`. To adjust, provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#'
#'  \item{car_rho}{The CAR model also has a spatial autocorrelation parameter, `rho`, which is assigned a uniform prior distribution. You can set the boundaries of the prior with: `ME$prior$car_rho <- uniform(lower_bound, upper_bound)`. You must specify values that are within the permissible range of values for `rho`; these are automatically printed to the console by \code{\link[geostan]{prep_car_data}}.}
#' }
#' }
#' }
#' 
#'
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. If `centerx` is a numeric vector, predictors will be centered on the values provided. This argument is passed to \code{\link[base]{scale}}.
#' 
#' @param prior_only Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.
#' @param chains Number of MCMC chains to use. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details. 
#' 
#' @param ... Other arguments passed to \code{\link[rstan]{sampling}}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' 
#' @details
#'
#' CAR models are discussed in Cressie and Wikle (2011, p. 184-88), Cressie (2015, Ch. 6-7), and Haining and Li (2020, p. 249-51). 
#'
#' The Stan code for this implementation of the CAR model first introduced in Donegan et al. (2021, supplementary material) for models of small area survey data.
#' 
#' Details and results depend on the \code{family} argument, as well as on the particular CAR specification chosen (see \link[geostan]{prep_car_data}).
#' 
#' ###  `auto_gaussian()`
#'
#' With \code{family = auto_gaussian()} (referred to here as the auto-Gaussian model), the CAR model is specified as follows:
#' ```
#'         Y ~ MVGauss(Mu, Sigma)
#'         Sigma = (I - rho C)^-1 * M * tau^2
#' ```
#' where \code{Mu} is the mean vector (with intercept, covariates, etc.), \code{C} is a spatial connectivity matrix, and \code{M} is a known diagonal matrix with diagonal entries proportional to the conditional variances. `C` and `M` are provided by \code{\link[geostan]{prep_car_data}}.
#'
#' The covariance matrix of the CAR model, \code{Sigma}, contains two parameters: \code{car_rho} (rho), which controls the degree of spatial autocorrelation, and the scale parameter, \code{car_scale} (tau). The range of permissible values for \code{rho} depends on the specification of \code{C} and \code{M}; for options, see \code{\link[geostan]{prep_car_data}} and Cressie and Wikle (2011, pp. 184-188).
#'
#' The auto-Gaussian model contains an implicit spatial trend (i.e., autocorrelation) component which is calculated as follows (Cressie 2015, p. 564):
#' ```
#'         trend = rho * C * (Y - Mu).
#' ```
#' This term can be extracted from a fitted auto-Gaussian model using the \code{\link[geostan]{spatial}} method.
#'
#' When applied to a fitted auto-Gaussian model, the \code{\link[geostan]{residuals.geostan_fit}} method returns `de-trended' residuals by default. That is,
#' ```
#'         residual = Y - Mu - trend.
#' ```
#' To obtain "raw" residuals (`Y - Mu`), use `residuals(fit, detrend = FALSE)`.
#' 
#' ### `poisson()`
#'
#' For \code{family = poisson()}, the model is specified as:
#'
#' ```
#'         Y ~ Poisson(exp(offset + lambda))
#'         lambda ~ MVGauss(Mu, Sigma)
#'         Sigma = (I - rho C)^-1 * M * tau^2
#' ```
#' These models are most often used to calculate small area incidence rates (mortality or disease incidence rates); the user provided offset should be, then, the natural logarithm of the denominator in the rates, e.g., log-population at risk.
#' 
#' For Poisson models, the \code{\link[geostan]{spatial}} method returns the parameter vector \code{phi}, which is the log-risk minus the intercept and any covariates:
#'  ```
#'         phi = lambda - Mu.
#' ```
#' This is the spatial autocorrelation component. This is equivalent to specifying the model as:
#' ```
#'         Y ~ Poisson(exp(offset + Mu + phi))
#'         phi ~ MVGauss(0, Sigma)
#'         Sigma = (I - rho C)^-1 * M * tau^2.
#' ```
#' 
#' In the Poisson CAR model, `phi` contains a latent spatial trend as well as additional variation around it. If you would like to extract the latent/implicit spatial trend from \code{phi}, you can do so by calculating (following Cressie 2015, p. 564):
#' ```
#'         trend = rho * C * phi.
#' ```
#' 
#' ### `binomial()`
#' 
#' For `family = binomial()`, the model is specified as:
#'``` 
#'         Y ~ Binomial(N, theta)
#'         logit(theta) ~ MVGauss(Mu, Sigma)
#'         Sigma = (I - rho C)^-1 * M * tau^2
#'```
#' where outcome data `Y` are counts, `N` is the number of trials, and `theta` is the 'success' rate. Note that the model formula should be structured as: `cbind(sucesses, failures) ~ x`, such that `trials = successes = failures`.
#' 
#' For fitted Binomial models, the \code{\link[geostan]{spatial}} method will return the parameter vector \code{phi}, equivalent to:
#'```
#'         phi = logit(theta) - Mu.
#'```
#'
#' ### Spatially lagged covariates (SLX)
#' 
#' The `slx` argument is a convenience function for including SLX terms. For example,
#' ```
#' stan_glm(y ~ x1 + x2, slx = ~ x1, \...)
#' ```
#' is a shortcut for
#' ```
#' stan_glm(y ~ I(W \%*\% x1) + x1 + x2, \...)
#' ```
#' where `W` is a row-standardized spatial weights matrix (see \code{\link[geostan]{shape2mat}}). SLX terms will always be *prepended* to the design matrix, as above, which is important to know when setting prior distributions for regression coefficients.
#'
#' For measurement error (ME) models, the SLX argument is the only way to include spatially lagged covariates since the SLX term needs to be re-calculated on each iteration of the MCMC algorithm.
#' 
#' ### Measurement error (ME) models
#' 
#' The ME models are designed for surveys with spatial sampling designs, such as the American Community Survey (ACS) estimates. With estimates, `x`, and their standard errors, `se`, the ME models have one of the the following two specifications, depending on the user input:
#' ```
#'        x ~ Gauss(x_true, se)
#'        x_true ~ MVGauss(mu, Sigma)
#'        Sigma = (I - rho C)^(-1) M * tau^2
#'        mu ~ Gauss(0, 100)
#'        tau ~ student_t(10, 0, 40)
#'        rho ~ uniform(lower_bound, upper_bound)
#' ```
#' where the covariance matrix, `Sigma`, has the conditional autoregressive specification, and `tau` is the scale parameter. If `ME$car_parts` is not provided by the user, then a non-spatial model will be used instead:
#' ```
#'        x ~ Gauss(x_true, se)
#'        x_true ~ student_t(df, mu, sigma)
#'        df ~ gamma(3, 0.2)
#'        mu ~ Gauss(0, 100)
#'        sigma ~ student_t(10, 0, 40)
#' ```
#' The observational error is the difference between the survey estimate and `x_true`, the actual value of the variable during the period of the survey.
#'
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame.}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including CAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the varying intercepts (\code{re}) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' 
#' \item{x_center}{If covariates are centered internally (i.e., `centerx` is not `FALSE`), then `x_centers` is the numeric vector of values on which the covariates were centered.}
#' \item{spatial}{A data frame with the name of the spatial component parameter (either "phi" or, for auto Gaussian models, "trend") and method ("CAR")}
#' \item{ME}{A list indicating if the object contains an ME model; if so, the user-provided ME list is also stored here.}
#' \item{C}{Spatial connectivity matrix (in sparse matrix format).}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source
#'
#' Cressie, Noel (2015 (1993)). *Statistics for Spatial Data*. Wiley Classics, Revised Edition.
#' 
#' Cressie, Noel and Wikle, Christopher (2011). *Statistics for Spatio-Temporal Data*. Wiley.
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Donegan, Connor (2021). Spatial conditional autoregressive models in Stan. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#' 
#' Haining, Robert and Li, Guangquan (2020). *Modelling Spatial and Spatial-Temporal Data: A Bayesian Approach*. CRC Press.
#' 
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(bayesplot)
#' options(mc.cores = parallel::detectCores())
#' data(sentencing)
#'
#' C <- shape2mat(sentencing, style = "B")
#' car.dl <- prep_car_data(C, style = "WCAR")
#' log_e <- log(sentencing$expected_sents)
#' fit.car <- stan_car(sents ~ offset(log_e),
#'                     family = poisson(),
#'                     data = sentencing,
#'                     car_parts = car.dl
#' )
#'
#' # MCMC diagnostics
#' rstan::stan_rhat(fit.car$stanfit)
#' rstan::stan_ess(fit.car$stanfit)
#'
#' # Spatial diagnostics
#' sp_diag(fit.car, sentencing)
#' 
#' # posterior predictive distribution
#' yrep <- posterior_predict(fit.car, S = 75)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep)
#'
#' # examine posterior distributions of CAR parameters
#' plot(fit.car, pars = c("car_scale", "car_rho"))
#'
#' # map the spatial autocorrelation term, phi
#' sp.trend <- spatial(fit.car)$mean
#' st_as_sf(sentencing) %>%
#'   ggplot() +
#'   geom_sf(aes(fill = sp.trend)) +
#'   scale_fill_gradient2()
#'  
#' # calculate log-standardized sentencing ratios (log-SSRs)
#' # (like Standardized Incidence Ratios: observed/exected case counts)
#' SSR <- fitted(fit.car)$mean
#' log.SSR <- log( SSR, base = 2 )
#'
#' st_as_sf(sentencing) %>%
#'  ggplot() +
#'  geom_sf(aes(fill = log.SSR)) +
#'  scale_fill_gradient2(
#'    midpoint = 0,
#'    name = NULL,
#'    breaks = seq(-3, 3, by = 0.5)
#'    ) +
#'  labs(title = "Log-Standardized Sentencing Ratios",
#'       subtitle = "log( Fitted/Expected ), base 2"
#'       ) +
#'  theme_void() +
#'  theme(
#'    legend.position = "bottom",
#'    legend.key.height = unit(0.35, "cm"),
#'    legend.key.width = unit(1.5, "cm")
#'  )
#'
#' ## DCAR specification (inverse-distance based)
#' A <- shape2mat(georgia, "B")
#' D <- sf::st_distance(sf::st_centroid(georgia))
#' A <- D * A
#' cp <- prep_car_data(A, "DCAR", k = 1)
#'
#' # view distance against inverse distance: (A+gamma)^(-k) / max(C)
#' plot(as.numeric(A),
#'     as.numeric(cp$C)
#'     )
#'
#' fit <- stan_car(log(rate.male) ~ 1,
#'                data = georgia,
#'                car = cp,
#'                cores = 5)
#'
#' me_diag(fit, georgia)
#' }
#'
#' @export
#' @md
#' @importFrom rstan extract_sparse_parts
#' 
stan_car <- function(formula,
                     slx,
                     re,
                     data,
                     car_parts,
                     family = gaussian(),
                     prior = NULL,                      
                     ME = NULL,                     
                     centerx = FALSE,
                     prior_only = FALSE,
                     chains = 4,
                     iter = 2e3,
                     refresh = 500,
                     pars = NULL,
                     control = NULL,
                     ...
                     ) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("gaussian", "auto_gaussian", "poisson", "binomial"))
    if (family$family == "gaussian") family <- auto_gaussian()
    stopifnot(!missing(data))
    check_car_parts(car_parts)
    stopifnot(length(car_parts$Delta_inv) == nrow(data))
    C <- car_parts$C        
    a.zero <- as.array(0, dim = 1)
    tmpdf <- as.data.frame(data)
    mod.mat <- model.matrix(formula, tmpdf)
    if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")  
    n <- nrow(mod.mat)
    family_int <- family_2_int(family)
    intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0) 
    if (intercept_only) {
        if (!missing(slx)) {
            stop("You provided a spatial lag of X (slx) term for an intercept only model. Did you intend to include a covariate?")
        }
        x_full <- x_no_Wx <- model.matrix(~ 0, data = tmpdf) 
        slx <- " "
        W.list <- list(w = 1, v = 1, u = 1)
        dwx <- 0
        wx_idx <- a.zero
  } else {
      xraw <- model.matrix(formula, data = tmpdf)
      xraw <- remove_intercept(xraw)
      x_no_Wx <- center_x(xraw, center = centerx)
      if (missing(slx)) {
          slx <- " "
          W.list <- list(w = 1, v = 1, u = 1)
          wx_idx = a.zero
          dwx <- 0
          x_full <- x_no_Wx          
      } else {
          stopifnot(inherits(slx, "formula"))
          W <- row_standardize(C, msg =  "Row standardizing connectivity matrix to calculate spatially lagged covaraite(s)")
          W.list <- rstan::extract_sparse_parts(W)
          Wx <- SLX(f = slx, DF = tmpdf, x = x_no_Wx, W = W)
          dwx <- ncol(Wx)
          wx_idx <- as.array( which(paste0("w.", colnames(x_no_Wx)) %in% colnames(Wx)), dim = dwx )
          x_full <- cbind(Wx, x_no_Wx)
      }
  }
    ModData <- make_data(formula, tmpdf, x_full)
    frame <- model.frame(formula, tmpdf)
    y <- y_int <- model.response(frame)
    if (family_int %in% c(1,2,5)) y_int <- rep(0, length(y))
    if (is.null(model.offset(frame))) {
        offset <- rep(0, times = n)
    } else {
        offset <- model.offset(frame)
    }
    if(missing(re)) {
        has_re <- n_ids <- id <- 0;
        id_index <- to_index(id, n = nrow(tmpdf))
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
        trials = rep(0, length(y)),
        n = n,
        offset = offset,
        has_re = has_re,
        n_ids = n_ids,
        id = id_index$idx,
        ## slx data -------------    
        W_w = as.array(W.list$w),
        W_v = as.array(W.list$v),
        W_u = as.array(W.list$u),
        dw_nonzero = length(W.list$w),
        dwx = dwx,
        wx_idx = wx_idx
    )
    ## PRIORS & CAR DATA -------------  
    is_student <- family$family == "student_t"
    priors_made <- make_priors(user_priors = prior,
                               y = y,
                               x = x_full,
                               link = family$link,
                               offset = offset)
    if (is.null(prior$car_scale)) {
        priors_made$car_scale <- priors_made$sigma        
    } else {
        priors_made$car_scale <- prior$car_scale
    }
    if (is.null(prior$car_rho)) {
        lims <- 1 / range(car_parts$lambda)
        priors_made$car_rho <- uniform(lims[1], lims[2])
    } else {
        priors_made$car_rho <- prior$car_rho
    }
    standata$car_rho_lims = c(priors_made$car_rho$lower, priors_made$car_rho$upper)     
    standata <- c(standata, car_parts)
    standata <- append_priors(standata, priors_made)    
    ## ME MODEL -------------
    # overwrites any user specified ME$car_parts
    if ( !is.null(ME) ) {
        if (is.logical(ME$car_parts) && !ME$car_parts) ME$car_parts <- NULL else ME$car_parts <- car_parts
    }    
    me.list <- prep_me_data(ME, x_no_Wx)
    standata <- c(standata, me.list)
    ## INTEGER OUTCOMES -------------    
    if (family$family == "binomial") {
        standata$y <- standata$y_int <- y[,1]
        standata$trials <- y[,1] + y[,2]
    }
    ## PARAMETERS TO KEEP, with CAR PARAMETERS [START] -------------            
    pars <- c(pars, 'intercept', 'car_scale', 'car_rho', 'fitted', 'log_lik')
    if (family_int < 5) pars <- c(pars, 'phi') 
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
    priors_made_slim <- priors_made[which(names(priors_made) %in% pars)]
    ## PARAMETERS TO KEEP, with CAR PARAMETERS [STOP] -------------                
    ## PRINT PRIORS -------------
    if (me.list$has_me) {
        if (me.list$spatial_me) {
            priors_made_slim <- c(priors_made_slim, list(ME_car_rho = me.list$ME_prior_car_rho))
        } else {
            priors_made_slim <- c(priors_made_slim, list(ME_df = me.list$ME_prior_df))
        }        
        priors_made_slim <- c(priors_made_slim, list(ME_location = me.list$ME_prior_mean, ME_scale = me.list$ME_prior_scale))
    }
    print_priors(prior, priors_made_slim)
    ## CALL STAN -------------  
    samples <- rstan::sampling(stanmodels$car, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
    ## OUTPUT -------------
    out <- clean_results(samples, pars, is_student, has_re, Wx, x_no_Wx, me.list$x_me_idx)
    out$data <- ModData
    out$family <- family
    out$formula <- formula
    out$slx <- slx
    out$re <- re_list
    out$priors <- priors_made_slim
    out$x_center <- attributes(x_full)$`scaled:center`
    out$ME <- list(has_me = me.list$has_me, spatial_me = me.list$spatial_me)
    if (out$ME$has_me) out$ME <- c(out$ME, ME)
    if (family_int == 5) {
        out$spatial <- data.frame(par = "trend", method = "CAR") #!#
    } else {
        out$spatial <- data.frame(par = "phi", method = "CAR")
    }
    out$C <- Matrix::Matrix(C)    
    R <- resid(out, summary = FALSE)
    out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C) )    
    return (out)
}

