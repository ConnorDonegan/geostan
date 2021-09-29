#' Intrinsic autoregressive models
#'
#' @description The intrinsic conditional auto-regressive (ICAR) model for spatial count data. Options include the BYM model, the BYM2 model, and a solo ICAR term. 
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' 
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms. 
#' 
#' @param re To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:
#' ```
#'        alpha_re ~ N(0, alpha_tau)
#'        alpha_tau ~ Student_t(d.f., location, scale).
#' ```
#' Before using this term, read the \code{Details} section and the \code{type} argument. Specifically, if you use `type = bym`, then an observational-level `re` term is already included in the model. (Similar for `type = bym2`.)
#' 
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'
#' @param C Spatial connectivity matrix which will be used to construct an edge list for the ICAR model, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms. It will automatically be row-standardized before calculating \code{slx} terms. \code{C} must be a binary symmetric \code{n x n} matrix.
#' 
#' @param type Defaults to "icar" (partial pooling of neighboring observations through parameter \code{phi}); specify "bym" to add a second parameter vector \code{theta} to perform partial pooling across all observations; specify "bym2" for the innovation introduced by Riebler et al. (2016). See \code{Details} for more information.
#' 
#' @param scale_factor For the BYM2 model, optional. If missing, this will be set to a vector of ones. See `Details`.
#'
#' @param family The likelihood function for the outcome variable. Current options are \code{binomial(link = "logit")} and \code{poisson(link = "log")}.
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
#' \item{nu}{`nu` is the degrees of freedom parameter in the Student's t likelihood (only used when `family = student_t()`). `nu` is assigned a gamma prior distribution. The default prior is `prior = list(nu = gamma(alpha = 3, beta = 0.2))`.
#' }
#'
#' \item{tau}{The scale parameter for random effects, or varying intercepts, terms. This scale parameter, `tau`, is assigned a half-Student's t prior. To set this, use, e.g., `prior = list(tau = student_t(df = 20, location = 0, scale = 20))`.}
#' }
#'
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list of priors. See the Details section below for more information. Elements of the \code{ME} list may include:
#' \describe{
#' 
#' \item{se}{A required dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' 
#' \item{bounds}{An optional numeric vector of length two providing the upper and lower bounds, respectively, of the variables. If not provided, they will be set to `c(-Inf, Inf)` (i.e., unbounded). Common usages include keeping percentages between zero and one hundred or proportions between zero and one.}
#'
#' \item{car_parts}{A list of data for the CAR model, as returned by \code{\link[geostan]{prep_car_data}}. If not provided, a non-spatial Student's t model will be used instead of the CAR model.}
#' 
#'  \item{prior}{Optionally provide a named list of prior distributions for the measurement error model(s) (see \code{\link[geostan]{priors}}). If none are provided, default priors will be assigned and printed to the console. \describe{
#' 
#'  \item{df}{If using a non-spatial ME model, the degrees of freedom for the Student's t model is assigned a gamma prior with default parameters of `gamma(alpha = 3, beta = 0.2)`. Provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#'
#'  \item{car_rho}{The CAR model also has a spatial autocorrelation parameter, `rho`, which is assigned a uniform prior distribution. You can set the boundaries of the prior with: `ME$prior$car_rho <- uniform(lower_bound, upper_bound)`. You must specify values that are within the permissible range of values for `rho`. This range is automatically printed to the console by \code{\link[geostan]{prep_car_data}}.}
#' 
#'   \item{location}{The prior for the location parameter, mu, is Gaussian (the default being `normal(location = 0, scale = 100)`). To adjust the prior distributions, provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#' 
#'  \item{scale}{The prior for the `scale` parameters is Student's t, and the default is `student_t(df = 10, location = 0, scale = 40)`. To adjust, provide values for each covariate in `se`, listing the values in the same order as the columns of `se`.}
#' }
#' }
#' }
#'
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. If `centerx` is a numeric vector, predictors will be centered on the values provided. This argument is passed to \code{\link[base]{scale}}.
#' 
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. 
#' @param iter Number of samples per chain. .
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details. 
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details
#' 
#'  The Stan code for the ICAR component of the model and the BYM2 option is from Morris et al. (2019) with adjustments to enable non-binary weights and disconnected graph structures (see Freni-Sterrantino (2018) and Donegan (2021)).
#'
#' The exact specification depends on the `type` argument. 
#'
#' ### 'icar'
#'
#' For Poisson models for count data, y, the basic model specification (`type = "icar"`) is:
#' ```
#'         y ~ Poisson(exp(offset + mu + phi))
#'         phi ~ ICAR(spatial_scale)
#'         spatial_scale ~ Gaussian(0, 1)
#' ```
#'  where `mu` contains an intercept and potentially covariates. The spatial trend, `phi`, has a mean of zero and a single scale parameter, `spatial_scale`.
#' 
#' The ICAR prior model is a CAR model that has a spatial autocorrelation parameter \code{car_alpha} equal to 1 (see \link[geostan]{stan_car}). Thus the ICAR prior places high probability on a smooth spatially (or temporally) varying mean. This is rarely sufficient to model the amount of variation present in social and health data.
#'
#' ### 'bym'
#'
#' Often, an observational-level random effect term, `theta`, is added to capture (heterogeneous or unstructured) deviations from `mu + phi`. The combined term is referred to as a convolution term:
#' ```
#'         convolution = phi + theta.
#' ```
#' This is known as the BYM model (Besag et al. 1991), and can be specified using `type = "bym"`:
#' ```
#'         y ~ Poisson(exp(offset + mu + phi + theta))
#'         phi ~ ICAR(spatial_scale)
#'         theta ~ Gaussian(0, theta_scale)
#'         spatial_scale ~ Gaussian(0, 1)
#'         theta_scale ~ Gaussian(0, 1)
#' ```
#'
#' ### 'bym2'
#' 
#' Riebler et al. (2016) introduce a variation on the BYM model (`type = "bym2"`). This specification combines `phi` and `theta` using a mixing parameter, `rho`, that controls the proportion of the variation that is attributable to the spatially autocorrelated term, `phi`, rather than the spatially unstructured term, `theta`. The terms share a single scale parameter:
#' ```
#'         convolution = [sqrt(rho/scale_factor) * phi_tilde + sqrt(1 - rho) * theta_tilde] * spatial_scale.
#'         phi_tilde ~ Gaussian(0, 1)
#'         theta_tilde ~ Gaussian(0, 1)
#'         spatial_scale ~ Gaussian(0, 1)
#' ```
#' The two `_tilde` terms are standard normal deviates, `rho` is restricted to values between zero and one, and `scale_factor` is a constant term provided by the user. By default, `scale_factor` is equal to one, so that it does nothing. Riebler et al. (2016) argue that the interpretation or meaning of the scale of the ICAR model depends on the graph structure, `C`. This implies that the same prior distribution assigned to the `spatial_scale` will differ in its implications if `C` is changed; in other words, the priors are not transportable across models, and models that use the same nominal prior actually have different priors assigned to `spatial_scale`.
#'
#' Borrowing `R` code from Morris (2017) and following Freni-Sterrantino et al. (2018), the following `R` code can be used to create the `scale_factor` for the BYM2 model (note, this requires the INLA R package), given a spatial adjacency matrix, C:
#' ```
#' ## create a list of data for stan_icar
#' icar.data <- geostan::prep_icar_data(C)
#' ## calculate scale_factor for each of k connected group of nodes
#' k <- icar.data$k
#' scale_factor <- vector(mode = "numeric", length = k)
#' for (j in 1:k) {
#'   g.idx <- which(icar.data$comp_id == j) 
#'   if (length(g.idx) == 1) {
#'        scale_factor[j] <- 1
#'        next
#'     }    
#'   Cg <- C[g.idx, g.idx] 
#'   scale_factor[j] <- scale_c(Cg) 
#' }
#' ```
#' This code adjusts for 'islands' or areas with zero neighbors, and it also handles disconnected graph structures (see Donegan 2021). Following Freni-Sterrantino (2018), disconnected components of the graph structure are given their own intercept term; however, this value is added to `phi` automatically inside the Stan model. Therefore, the user never needs to make any adjustments for this term. (If you want to avoid complications from a disconnected graph structure, see \code{\link[geostan]{stan_car}}).
#' 
#' Note, the code above requires the `scale_c` function; it has package dependencies that are not included in `geostan`. To use `scale_c`, you have to load the following `R` function:
#' ```
#' #' compute scaling factor for adjacency matrix, accounting for differences in spatial connectivity 
#' #'
#' #' @param C connectivity matrix
#' #'
#' #' @details
#' #'
#' #' Requires the following packages: 
#' #'
#' #' library(Matrix)
#' #' library(INLA);
#' #' library(spdep)
#' #' library(igraph)
#' #'  
#' #' @source
#' #'
#' #'   Morris, Mitzi (2017). Spatial Models in Stan: Intrinsic Auto-Regressive Models for Areal Data. <https://mc-stan.org/users/documentation/case-studies/icar_stan.html>
#' #'
#' scale_c <- function(C) {
#'  geometric_mean <- function(x) exp(mean(log(x))) 
#'  N = dim(C)[1]
#'  Q =  Diagonal(N, rowSums(C)) - C
#'  Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
#'  Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))
#'  scaling_factor <- geometric_mean(Matrix::diag(Q_inv)) 
#'  return(scaling_factor) 
#'}
#' ```
#'
#'  ### Spatially lagged covariates (SLX)
#' 
#' The `slx` argument is a convenience function for including SLX terms. For example,
#' ```
#' stan_glm(y ~ x1 + x2, slx = ~ x1, ...)
#' ```
#' is a shortcut for
#' ```
#' stan_glm(y ~ I(W %*% x1) + x1 + x2, ...)
#' ```
#' where `W` is a row-standardized spatial weights matrix (see \code{\link[geostan]{shape2mat}}). SLX terms will always be *prepended* to the design matrix, as above, which is important to know when setting prior distributions for regression coefficients.
#'
#' For measurement error (ME) models, the SLX argument is the only way to include spatially lagged covariates since the SLX term needs to be re-calculated on each iteration of the MCMC algorithm.
#' 
#' ### Measurement error (ME) models
#' 
#' The ME models are designed for surveys with spatial sampling designs, such as the American Community Survey (ACS) estimates (Donegan et al. 2021; Donegan 2021). With estimates, `x`, and their standard errors, `se`, the ME models have one of the the following two specifications, depending on the user input:
#' ```
#'        x ~ Gauss(x_true, se)
#'        x_true ~ MVGauss(mu, Sigma)
#'        Sigma = (I - rho * C)^(-1) M * tau^2
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
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors and the weight attached to each pair (i.e., their corresponding element in the connectivity matrix  C}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ICAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list with two name elements, \code{formula} and \code{Data}, containing the formula \code{re} and a data frame with columns \code{id} (the grouping variable) and \code{idx} (the index values assigned to each group).}
#' \item{priors}{Prior specifications.}
#' \item{x_center}{If covariates are centered internally (i.e., `centerx` is not `FALSE`), then `x_centers` is the numeric vector of values on which the covariates were centered.}
#' \item{spatial}{A data frame with the name of the spatial parameter (\code{"phi"} if \code{type = "icar"} else \code{"convolution"}) and method (\code{toupper(type)}).}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#'
#' @seealso \link[geostan]{shape2mat}, \link[geostan]{stan_car}, \link[geostan]{stan_esf}, \link[geostan]{stan_glm}, \link[geostan]{prep_icar_data}
#' 
#' @source
#'
#' Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society: Series B (Methodological), 36(2), 192-225.
#'
#' Besag, J., York, J., & Mollié, A. (1991). Bayesian image restoration, with two applications in spatial statistics. Annals of the institute of statistical mathematics, 43(1), 1-20.
#'
#' Donegan, Connor. 2021. Flexible functions for ICAR, BYM, and BYM2 models in Stan. Code repository. <https://github.com/ConnorDonegan/Stan-IAR>
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Donegan, Connor (2021). Spatial conditional autoregressive models in Stan. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#' 
#' Freni-Sterrantino, Anna, Massimo Ventrucci, and Håvard Rue. 2018. A Note on Intrinsic Conditional Autoregressive Models for Disconnected Graphs. Spatial and Spatio-Temporal Epidemiology 26: 25–34.
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. Spatial and spatio-temporal epidemiology, 31, 100301.
#'
#' Riebler, A., Sorbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Statistical Methods in Medical Research, 25(4), 1145-1165.
#'
#' @examples
#' \dontrun{
#' library(rstan)
#' library(bayesplot)
#' library(sf)
#' options(mc.cores = parallel::detectCores())
#' data(sentencing)
#'
#' C <- shape2mat(sentencing, "B")
#' log_e <- log(sentencing$expected_sents)
#' fit.bym <- stan_icar(sents ~ offset(log_e),
#'                      family = poisson(),
#'                      data = sentencing,
#'                      type = "bym",
#'                      C = C
#'  )
#'
#' # check effective sample size and convergence
#' rstan::stan_ess(fit.bym$stanfit)
#' rstan::stan_rhat(fit.bym$stanfit)
#'
#' # see some spatial diagnostics
#' sp_diag(fit.bym, sentencing)
#'
#' # posterior predictive distribution
#' yrep <- posterior_predict(fit.bym, S = 100)
#' y <- sentencing$sents
#' bayesplot::ppc_dens_overlay(y, yrep)
#' 
#' # map the smooth spatial term
#' sp.trend <- spatial(fit.bym)$mean
#' ggplot( st_as_sf(sentencing) ) +
#'   geom_sf(aes(fill = sp.trend)) +
#'   scale_fill_gradient2() +
#'   theme_void()
#'
#' # calculate log-standardized sentencing ratios (log-SSRs)
#' ## (like Standardized Incidence Ratios: observed/exected case counts)
#' f <- fitted(fit.bym)$mean
#' SSR <- f / sentencing$expected_sents
#' log.SSR <- log( SSR, base = 2)
#'
#' ggplot( st_as_sf(sentencing) ) +
#'   geom_sf(aes(fill = log.SSR)) +
#'   scale_fill_gradient2() +
#'   labs(title = "Log-standardized sentencing ratios",
#'        subtitle = "log( Fitted/Expected), base 2") +
#'   theme_void() +
#'   theme(
#'    legend.position = "bottom",
#'    legend.key.height = unit(0.35, "cm"),
#'    legend.key.width = unit(1.5, "cm")
#'   )
#' }
#'
#' @export
#' @md
#' @importFrom rstan extract_sparse_parts
stan_icar <- function(formula,
                      slx,
                      re,
                      data,
                      C,
                      family = poisson(),                      
                      type = c("icar", "bym", "bym2"),                      
                      scale_factor = NULL,
                      prior = NULL,                      
                      ME = NULL,
                      centerx = FALSE,
                      prior_only = FALSE,
                      chains = 4, iter = 2e3, refresh = 500,
                      pars = NULL,
                      control = NULL,
                      ...) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("poisson", "binomial"))
    stopifnot(!missing(data))
    stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
    stopifnot(all(dim(C) == nrow(data)))
    #### ICAR TYPE [START] --------
    type <- match.arg(type)
    #### ICAR TYPE [STOP] --------
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
    if (family_int %in% c(1,2)) y_int <- rep(0, length(y))
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
    ## PRIORS -------------  
    is_student <- family$family == "student_t"
    priors_made <- make_priors(user_priors = prior,
                               y = y,
                               x = x_full,
                               link = family$link,
                               offset = offset)
    standata <- append_priors(standata, priors_made)        
    ## ICAR DATA [START] -------------
    iar.list <- prep_icar_data(C, scale_factor = scale_factor)
    if (!inherits(scale_factor, "NULL")) stopifnot(length(scale_factor) == iar.list$k)
    standata <- c(standata, iar.list)
    standata$type <- match(type, c("icar", "bym", "bym2"))
    ## ICAR DATA [STOP] -------------    
    ## EMPTY PLACEHOLDERS
    standata <- c(standata, empty_esf_data(n))
    ## ME MODEL -------------  
    me.list <- prep_me_data(ME, x_no_Wx)
    standata <- c(standata, me.list)  
    ## INTEGER OUTCOMES -------------    
    if (family$family == "binomial") {
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
    }
    ## PARAMETERS TO KEEP with ICAR [START] -------------        
    pars <- c(pars, 'intercept', 'log_lik', 'fitted', 'phi', 'spatial_scale')
    if (type == "bym2") pars <- c(pars, "theta", "rho")
    if (type == "bym") pars <- c(pars, "theta", "theta_scale")
    if (standata$m) pars <- c(pars, "alpha_phi")
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
    if (me.list$has_me) {
        if (me.list$spatial_me) {
            priors_made_slim <- c(priors_made_slim, list(ME_car_rho = me.list$ME_prior_car_rho))
        } else {
            priors_made_slim <- c(priors_made_slim, list(ME_df = me.list$ME_prior_df))
        }
        priors_made_slim <- c(priors_made_slim, list(ME_location = me.list$ME_prior_mean, ME_scale = me.list$ME_prior_scale))
    }
    print_priors(prior, priors_made_slim)    
    ## PARAMETERS TO KEEP with ICAR [STOP] -------------              
    ## CALL STAN -------------  
    samples <- rstan::sampling(stanmodels$foundation, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
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
    ## ICAR OUTPUT [START] --------
    out$edges <- edges(C)       
    out$spatial <- data.frame(par = "phi", method = toupper(type))
    ## ICAR OUTPUT [STOP] --------
    R <- resid(out, summary = FALSE)
    out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C) )    
    return (out)
}

