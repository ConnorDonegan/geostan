#' Conditional autoregressive (CAR) models
#' 
#' @description Use the CAR model as a prior on parameters, or fit data to a spatial Gaussian CAR model.
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
#' With the CAR model, any \code{alpha_re} term should be at a *different* level or scale than the observations; that is, at a different scale than the autocorrelation structure of the CAR model itself.
#'
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' 
#' @param car_parts A list of data for the CAR model, as returned by \code{\link[geostan]{prep_car_data}}.
#'
#' @param C Optional spatial connectivity matrix which will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms; it will automatically be row-standardized before calculating \code{slx} terms. See \code{\link[geostan]{shape2mat}}.
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
#' \item{car_scale}{Scale parameter for the CAR model, \code{car_scale}. The scale is assigned a Student's t prior model (constrained to be positive).}
#'
#' \item{car_rho}{The spatial autocorrelation parameter in the CAR model, `rho`, is assigned a uniform prior distribution. By default, the prior will be uniform over all permissible values as determined by the eigenvalues of the connectivity matrix, `C`. The range of permissible values for `rho` is automatically printed to the console by \code{\link[geostan]{prep_car_data}}.}
#'
#' \item{tau}{The scale parameter for any varying intercepts (a.k.a exchangeable random effects, or partial pooling) terms. This scale parameter, `tau`, is assigned a Student's t prior (constrained to be positive).}
#' 
#' }
#' 
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a list of data as constructed by the \code{\link[geostan]{prep_me_data}} function. 
#'
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. If the ME argument is used, the modeled covariate (i.e., latent variable), rather than the raw observations, will be centered. When using the ME argument, this is the recommended method for centering the covariates.
#'
#' @param censor_point Integer value indicating the maximum censored value; this argument is for modeling censored (suppressed) outcome data, typically disease case counts or deaths. 
#' 
#' @param prior_only Logical value; if \code{TRUE}, draw samples only from the prior distributions of parameters.
#' @param chains Number of MCMC chains to use. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param keep_all  If `keep_all = TRUE` then samples for all parameters in the Stan model will be kept; this is necessary if you want to do model comparison with Bayes factors and the `bridgesampling` package.
#' @param control A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details. 
#' 
#' @param ... Other arguments passed to \code{\link[rstan]{sampling}}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' 
#' @details
#'
#' CAR models are discussed in Cressie and Wikle (2011, p. 184-88), Cressie (2015, Ch. 6-7), and Haining and Li (2020, p. 249-51). It is often used for areal or lattice data.
#'
#' Details for the Stan code for this implementation of the CAR model can be found in Donegan (2021).
#'
#' The general scheme for the CAR model is as follows:
#' \deqn{
#'  y \sim Gauss( \mu, ( I - \rho C)^{-1} M),
#' }
#' where \eqn{I} is the identity matrix, \eqn{\rho} is a spatial dependence parameter, \eqn{C} is a spatial connectivity matrix, and \eqn{M} is a diagonal matrix of variance terms. The diagonal of \eqn{M} contains a scale parameter \eqn{\tau} multiplied by a vector of weights (often set to be proportional to the inverse of the number of neighbors assigned to each site). The CAR model owes its name to the fact that this joint distribution corresponds to a set of conditional distributions that relate the expected value of each observation to a function of neighboring values, i.e., the Markov condition holds:
#' \deqn{
#' E(y_i | y_1, y_2, \dots y_{i-1}, y_{i+1}, y_n) = \mu_i + \rho \sum_{j=1}^n c_{i,j} (y_j - \mu_j),
#' }
#' where entries of \eqn{c_{i,j}} are non-zero only if \eqn{j \in N(i)}, and \eqn{N(i)} indexes the sites that are neighbors of the \eqn{i^{th}} site.
#' 
#' With the Gaussian probability distribution,
#' \deqn{
#'  y_i | y_j: j \neq i \sim Gauss(\mu_i + \rho \sum_{j=1}^n c_{i,j} (y_j - \mu_j), \tau_i^2)
#' }
#' where \eqn{\tau_i} is a scale parameter and \eqn{\mu_i} may contain covariates or simply the intercept.
#'
#' The covariance matrix of the CAR model contains two parameters: \code{car_rho} (\eqn{\rho}), which controls the degree of spatial autocorrelation, and the scale parameter, \code{car_scale} (\eqn{\tau}). The range of permissible values for \eqn{\rho} depends on the specification of \eqn{\boldsymbol C} and \eqn{\boldsymbol M}; for specification options, see \code{\link[geostan]{prep_car_data}} and Cressie and Wikle (2011, pp. 184-188) or Donegan (2021).
#' 
#' Further details of the models and results depend on the \code{family} argument, as well as on the particular CAR specification chosen (see \link[geostan]{prep_car_data}).
#'
#' ###  Auto-Gaussian
#'
#' When \code{family = auto_gaussian()} (the default), the CAR model is applied directly to the data as follows:
#' \deqn{
#'  y \sim Gauss( \mu, (I - \rho C)^{-1} M),
#' }
#' where \eqn{\mu} is the mean vector (with intercept, covariates, etc.), \eqn{C} is a spatial connectivity matrix, and \eqn{M} is a known diagonal matrix containing the conditional variances \eqn{\tau_i^2}. `C` and `M` are provided by \code{\link[geostan]{prep_car_data}}.
#'
#' The auto-Gaussian model contains an implicit spatial trend (i.e., autocorrelation) component \eqn{\phi} which can be calculated as follows (Cressie 2015, p. 564):
#' \deqn{
#'  \phi = \rho C (y - \mu).
#' }
#' This term can be extracted from a fitted auto-Gaussian model using the \code{\link[geostan]{spatial}} method.
#'
#' When applied to a fitted auto-Gaussian model, the \code{\link[geostan]{residuals.geostan_fit}} method returns 'de-trended' residuals \eqn{R} by default. That is,
#' \deqn{
#' R = y - \mu - \rho C (y - \mu).
#' }
#' To obtain "raw" residuals (\eqn{y - \mu}), use `residuals(fit, detrend = FALSE)`.
#' 
#' ### Poisson
#'
#' For \code{family = poisson()}, the model is specified as:
#'\deqn{
#' y \sim Poisson(e^{O + \lambda}) \\
#' \lambda \sim Gauss(\mu, (I - \rho C)^{-1} \boldsymbol M).
#' }
#' If the raw outcome consists of a rate \eqn{\frac{y}{p}} with observed counts \eqn{y} and denominator {p} (often this will be the size of the population at risk), then \eqn{O=log(p)} is the log of the denominator (the offset term).
#'
#' This is often written (equivalently) as:
#' \deqn{
#' y \sim Poisson(e^{O + \mu + \phi}) \\
#' \phi \sim Gauss(0, (I - \rho C)^{-1} \boldsymbol M).
#' }
#' For Poisson models, the \code{\link[geostan]{spatial}} method returns the parameter vector \eqn{\phi}.
#' 
#' In the Poisson CAR model, \eqn{\phi} contains a latent spatial trend as well as additional variation around it: \eqn{\phi_i = \rho \sum_{i=1}^n c_{ij} \phi_j + \epsilon_i}, \eqn{\epsilon_i \sim Gauss(0, \tau_i^2)}. If you would like to extract the latent/implicit spatial trend from \eqn{\phi}, you can do so by calculating (following Cressie 2015, p. 564):
#' \deqn{
#' \rho  C  \phi.
#' }
#' 
#' ### Binomial
#' 
#' For `family = binomial()`, the model is specified as:
#'\deqn{
#' y \sim Binomial(N, \lambda)  \\
#' logit(\lambda) \sim Gauss(\mu, (I - \rho C)^{-1} \boldsymbol M).
#' }
#' where outcome data `y` are counts, `N` is the number of trials, and `theta` is the 'success' rate. Note that the model formula should be structured as: `cbind(sucesses, failures) ~ x`, such that `trials = successes + failures`.
#'
#' This is often written (equivalently) as:
#' \deqn{
#' y \sim Binomial(N, (\mu + \phi))  \\
#' logit(\phi) \sim Gauss(0, (I - \rho C)^{-1} \boldsymbol M).
#' }
#' For fitted Binomial models, the \code{\link[geostan]{spatial}} method will return the parameter vector \code{phi}.
#' 
#' As is also the case for the Poisson model, \eqn{\phi} contains a latent spatial trend as well as additional variation around it. If you would like to extract the latent/implicit spatial trend from \eqn{\phi}, you can do so by calculating:
#' \deqn{
#' \rho C \phi.
#' }
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
#' The ME models are designed for surveys with spatial sampling designs, such as the American Community Survey (ACS) estimates. Given estimates \eqn{x}, their standard errors \eqn{s}, and the target quantity of interest (i.e., the unknown true value) \eqn{z}, the ME models have one of the the following two specifications, depending on the user input. If a spatial CAR model is specified, then:
#' \deqn{
#'  x \sim Gauss(z, s^2) \\
#'  z \sim Gauss(\mu_z, \Sigma_z) \\
#' \Sigma_z = (I - \rho C)^{-1} M \\
#'  \mu_z \sim Gauss(0, 100) \\
#'  \tau_z \sim Student(10, 0, 40), \tau > 0 \\
#'  \rho_z \sim uniform(l, u)
#'  }
#' where \eqn{\Sigma} specifies a spatial conditional autoregressive model with scale parameter \eqn{\tau} (on the diagonal of \eqn{M}), and \eqn{l}, \eqn{u} are the lower and upper bounds that \eqn{\rho} is permitted to take (which is determined by the extreme eigenvalues of the spatial connectivity matrix \eqn{C}).
#' 
#' For non-spatial ME models, the following is used instead:
#' \deqn{
#' x \sim Gauss(z, s^2) \\
#' z \sim student(\nu_z, \mu_z, \sigma_z) \\
#' \nu_z \sim gamma(3, 0.2) \\
#' \mu_z \sim Gauss(0, 100) \\
#' \sigma_z \sim student(10, 0, 40).
#' }
#' 
#' For strongly skewed variables, such as census tract poverty rates, it can be advantageous to apply a logit transformation to \eqn{z} before applying the CAR or Student-t prior model. When the `logit` argument is used, the model becomes:
#' \deqn{
#' x \sim Gauss(z, s^2) \\
#' logit(z) \sim Gauss(\mu_z, \Sigma_z) 
#' ...
#' }
#' and similarly for the Student t model:
#' \deqn{
#' x \sim Gauss(z, s^2) \\
#' logit(z) \sim student(\nu_z, \mu_z, \sigma_z) \\
#' ...
#' }
#'
#' ### Censored counts
#'
#' Vital statistics systems and disease surveillance programs typically suppress case counts when they are smaller than a specific threshold value. In such cases, the observation of a censored count is not the same as a missing value; instead, you are informed that the value is an integer somewhere between zero and the threshold value. For Poisson models (`family = poisson())`), you can use the `censor_point` argument to encode this information into your model. 
#'
#' Internally, `geostan` will keep the index values of each censored observation, and the index value of each of the fully observed outcome values. For all observed counts, the likelihood statement will be:
#' \deqn{
#' p(y_i | \text{data}, \text{model}) = poisson(y_i | \mu_i), 
#' }
#' as usual, where \eqn{\mu_i} may include whatever spatial terms are present in the model.
#'
#' For each censored count, the likelihood statement will equal the cumulative Poisson distribution function for values zero through the censor point:
#' \deqn{
#' p(y_i | \text{data}, \text{model}) = \sum_{m=0}^{M} Poisson( m | \mu_i),
#' }
#' where \eqn{M} is the censor point and \eqn{\mu_i} again is the fitted value for the \eqn{i^{th}} observation.
#' 
#' For example, the US Centers for Disease Control and Prevention's CDC WONDER database censors all death counts between 0 and 9. To model CDC WONDER mortality data, you could provide `censor_point = 9` and then the likelihood statement for censored counts would equal the summation of the Poisson probability mass function over each integer ranging from zero through 9 (inclusive), conditional on the fitted values (i.e., all model parameters). See Donegan (2021) for additional discussion, references, and Stan code.
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
#' \item{x_center}{If covariates are centered internally (`centerx = TRUE`), then `x_center` is a numeric vector of the values on which covariates were centered.}
#' 
#' \item{spatial}{A data frame with the name of the spatial component parameter (either "phi" or, for auto Gaussian models, "trend") and method ("CAR")}
#' \item{ME}{A list indicating if the object contains an ME model; if so, the user-provided ME list is also stored here.}
#' \item{C}{Spatial connectivity matrix (in sparse matrix format).}
#' }
#' 
#' @author Connor Donegan, \email{connor.donegan@gmail.com}
#' 
#' @source
#'
#' Besag, Julian (1974). Spatial interaction and the statistical analysis of lattice systems. *Journal of the Royal Statistical Society* B36.2: 192–225.
#' 
#' Cressie, Noel (2015 (1993)). *Statistics for Spatial Data*. Wiley Classics, Revised Edition.
#' 
#' Cressie, Noel and Wikle, Christopher (2011). *Statistics for Spatio-Temporal Data*. Wiley.
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Donegan, Connor (2021). Building spatial conditional autoregressive (CAR) models in the Stan programming language. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#' 
#' Haining, Robert and Li, Guangquan (2020). *Modelling Spatial and Spatial-Temporal Data: A Bayesian Approach*. CRC Press.
#' 
#'
#' @examples
#' 
#' # model mortality risk
#' data(georgia)
#' C <- shape2mat(georgia, style = "B")
#' cp <- prep_car_data(C)
#' 
#' fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)),
#'                 car_parts = cp,
#'                 data = georgia,
#'                 family = poisson(),
#'                 iter = 800, chains = 2 # for example speed only
#'                  )
#' rstan::stan_rhat(fit$stanfit)
#' rstan::stan_mcse(fit$stanfit)
#' print(fit)
#' sp_diag(fit, georgia)
#'
#' # censored count outcomes
#' sum(is.na(georgia$deaths.female))
#' fit <- stan_car(deaths.female ~ offset(log(pop.at.risk.female)),
#'                 car_parts = cp,
#'                 data = georgia,
#'                 family = poisson(),
#'                 censor_point = 9,
#'                 iter = 800, chains = 2 # for example speed only
#'                 )
#' 
#' ## DCAR specification (inverse-distance based)
#' library(sf)
#' A <- shape2mat(georgia, "B")
#' D <- sf::st_distance(sf::st_centroid(georgia))
#' A <- D * A
#' cp <- prep_car_data(A, "DCAR", k = 1)
#' 
#' fit <- stan_car(deaths.male ~ offset(log(pop.at.risk.male)),
#'                data = georgia,
#'                car = cp,
#'                family = poisson(),
#'                iter = 800, chains = 1 # for example speed only 
#' )
#' print(fit)
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
                     C,
                     family = gaussian(),
                     prior = NULL,                      
                     ME = NULL,                     
                     centerx = FALSE,
                     prior_only = FALSE,
                     censor_point,                     
                     chains = 4,
                     iter = 2e3,
                     refresh = 500,
                     keep_all = FALSE,
                     pars = NULL,
                     control = NULL,
                     ...
                     ) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("gaussian", "auto_gaussian", "poisson", "binomial"))
    if (family$family == "gaussian" | family$family == "auto_gaussian") family <- auto_gaussian(type = "CAR")
    stopifnot(!missing(data))
    check_car_parts(car_parts)
    stopifnot(length(car_parts$Delta_inv) == nrow(data))
    if (!missing(C)) {
        stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
        stopifnot(all(dim(C) == nrow(data)))
    } else {
        C <- car_parts$C
        if (car_parts$WCAR == 0) {
            message("Consider providing the matrix C explicitly using the C argument. The matrix C is used for calculating spatial-lag of X (SLX) terms and residual spatial autocorrelation. Since you did not provide C, the matrix is being taken from car_parts$C.")
        }
    }
    tmpdf <- as.data.frame(data)
    n <- nrow(tmpdf)    
    family_int <- family_2_int(family)        
    if (!missing(censor_point)) if (family$family != "poisson") stop("censor_point argument is only available for Poisson models.")
    if (missing(censor_point)) censor_point <- FALSE
    mod_frame <- model.frame(formula, tmpdf, na.action = NULL)
    handle_missing_x(mod_frame)
    y_index_list <- handle_censored_y(censor_point, mod_frame)
    y <- y_int <- model.response(mod_frame)
    ## CAR: INCLUDE AUTO-GAUSSIAN AS family_int=5 -------------      
    if (family_int %in% c(1,2,5)) y_int <- rep(0, length(y))
    ## -------------      
    y[y_index_list$y_mis_idx] <- y_int[y_index_list$y_mis_idx] <- 0
    mod_frame[y_index_list$y_mis_idx, 1] <- 0 
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
          W <- row_standardize(C, msg =  "Row standardizing connectivity matrix to calculate spatially lagged covaraite(s)")
          W.list <- rstan::extract_sparse_parts(W)
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
        trials = rep(0, length(y)),
        n = n,
        offset = offset,
        has_re = has_re,
        n_ids = n_ids,
        id = id_index$idx,
        center_x = centerx,  ####////!!!!####        
        ## slx data -------------    
        W_w = as.array(W.list$w),
        W_v = as.array(W.list$v),
        W_u = as.array(W.list$u),
        dw_nonzero = length(W.list$w),
        dwx = dwx,
        wx_idx = wx_idx
    )
    ## ADD MISSING/OBSERVED INDICES -------------  
    standata <- c(y_index_list, standata)    
    ## PRIORS & CAR DATA -------------  
    is_student <- FALSE ##family$family == "student_t"
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
    standata$car <- 1
    ## EMPTY PLACEHOLDERS
    standata <- c(standata, empty_icar_data(n), empty_esf_data(n), empty_sar_data(n))    
    ## ME MODEL -------------
    me.list <- make_me_data(ME, xraw)
    standata <- c(standata, me.list)
    ## INTEGER OUTCOMES -------------    
    if (family$family == "binomial") {
        standata$y <- standata$y_int <- y[,1]
        standata$trials <- y[,1] + y[,2]
    }
    ## PARAMETERS TO KEEP, with CAR PARAMETERS [START] -------------            
    pars <- c(pars, 'intercept', 'car_scale', 'car_rho', 'fitted', 'log_lik')
    if (family_int < 5) pars <- c(pars, 'log_lambda_mu') 
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
    if (me.list$has_me) priors_made_slim$ME_model <- ME$prior    
    ## PRINT PRIORS  -------------
    print_priors(prior, priors_made_slim)
    ## MCMC INITIAL VALUES -------------
    inits <- "random"
    if (censor_point > 0 | (standata$has_me == 1 && any(standata$use_logit == 1))) {
        FRAME <- sys.nframe()
        inits <- init_fn_builder(FRAME_NUMBER = FRAME)
    }     
    ## CALL STAN -------------  
    standata$car <- 1
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
    if (family_int == 5) {
        out$spatial <- data.frame(par = "trend", method = "CAR") #!#
    } else {
        out$spatial <- data.frame(par = "phi", method = "CAR")
    }
    out$C <- Matrix::Matrix(C)    
    R <- resid(out, summary = FALSE)
    out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C, warn = FALSE, na.rm = TRUE) )    
    return (out)
}

