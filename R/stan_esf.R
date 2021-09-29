#' Spatial filtering
#'
#' @description Fit a spatial regression model using eigenvector spatial filtering (ESF).
#' 
#' @param formula A model formula, following the R \code{\link[stats]{formula}} syntax. Binomial models are specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#'
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms. 
#' 
#' @param re To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:
#' ```
#'        alpha_re ~ N(0, alpha_tau)
#'        alpha_tau ~ Student_t(d.f., location, scale).
#' ```
#' 
#' @param C Spatial connectivity matrix which will be used to calculate eigenvectors, if `EV` is not provided by the user. Typically, the binary connectivity matrix is best for calculating eigenvectors (i.e., using `C = shape2mat(shape, style = "B")`). This matrix will also be used to calculate residual spatial autocorrelation and any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms. See \code{\link[geostan]{shape2mat}}.
#'
#' @param nsa Include eigenvectors representing negative spatial autocorrelation? Defaults to \code{nsa = FALSE}. This is ignored if \code{EV} is provided.
#' 
#' @param threshold Eigenvectors with standardized Moran coefficient values below this `threshold` value will be excluded from the candidate set of eigenvectors, `EV`. This defaults to \code{threshold = 0.25}, and is ignored if \code{EV} is provided. 
#' 
#' @param EV A matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial (see \code{\link[geostan]{make_EV}}). If `EV` is provided, still also provide a spatial weights matrix \code{C} for other purposes; `threshold` and `nsa` are ignored for user provided `EV`.
#' 
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'
#' @param family The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{student_t()} and \code{poisson(link = "log")}, and \code{binomial(link = "logit")}.
#'
#' @param prior A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):  \describe{
#' 
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model `formula`. Note that if you also use `slx` terms (spatially lagged covariates), and you use custom priors for `beta`, then you have to provide priors for the slx terms. Since slx terms are *prepended* to the design matrix, the prior for the slx term will be listed first.
#' }
#'
#' \item{sigma}{For `family = gaussian()` and `family = student_t()` models, the scale parameter, `sigma`, is assigned a (half-) Student's t prior distribution. The half-Student's t prior for `sigma` is constrained to be positive.}
#'
#' \item{nu}{`nu` is the degrees of freedom parameter in the Student's t likelihood (only used when `family = student_t()`). `nu` is assigned a gamma prior distribution. The default prior is `prior = list(nu = gamma(alpha = 3, beta = 0.2))`. }
#'
#' \item{tau}{The scale parameter for random effects, or varying intercepts, terms. This scale parameter, `tau`, is assigned a half-Student's t prior. To set this, use, e.g., `prior = list(tau = student_t(df = 20, location = 0, scale = 20))`.}
#'
#' \item{beta_ev}{The eigenvector coefficients are assigned the horseshoe prior (Piironen and Vehtari, 2017), parameterized by `global_scale` (to control overall prior sparsity), plus the degrees of freedom and scale of a Student's t model for any large coefficients (see \code{\link[geostan]{priors}}). To allow the spatial filter to account for a greater amount of spatial autocorrelation (i.e., if you find the residuals contain spatial autocorrelation), increase the global scale parameter (to a maximum of `global_scale = 1`).}
#' }
#' 
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
#' 
#' @param centerx To center predictors on their mean values, use `centerx = TRUE`. If `centerx` is a numeric vector, predictors will be centered on the values provided. This argument is passed to \code{\link[base]{scale}}.
#' 
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. 
#' 
#' @param ... Other arguments passed to \link[rstan]{sampling}. 
#' @details
#'
#' Eigenvector spatial filtering (ESF) is a method for spatial regression analysis. ESF is extensively covered in Griffith et al. (2019). This function implements the methodology introduced in Donegan et al. (2020), which uses Piironen and Vehtari's (2017) regularized horseshoe prior.
#'
#' ESF decomposes spatial autocorrelation into a linear combination of various patterns, typically at different scales (such as local, regional, and global trends). By adding a spatial filter to a regression model, any spatial autocorrelation is shifted from the residuals to the spatial filter. ESF models take the spectral decomposition of a transformed spatial connectivity matrix, \code{C}. The resulting eigenvectors, `EV`, are mutually orthogonal and uncorrelated map patterns. The spatial filter is `EV * beta_ev`, where `beta_ev` is a vector of coefficients.
#'
#' ESF decomposes the data into a global mean, `alpha`, global patterns contributed by covariates, `X * beta`, spatial trends, `EV * beta_ev`, and residual variation. Thus, for `family=gaussian()`,
#' 
#' ```
#'         Y ~ Gauss(alpha + X * beta + EV * beta_ev, sigma).
#'```
#' An ESF component can be incorporated into the linear predictor of any generalized linear model. For example, a spatial Poisson model for rare disease incidence may be specified as follows:
#' ```
#'         Y ~ Poisson(exp(offset + Mu))
#'         Mu = alpha + EV * beta_ev + A
#'         A ~ Guass(0, tau)
#'         tau ~ student(20, 0, 2)
#'         beta_ev ~ horseshoe(.)
#' ```
#' 
#' The \code{\link[geostan]{spatial.geostan_fit}} method will return `EV * beta`.
#'
#' The model can also be extended to the space-time domain; see \link[geostan]{shape2mat} to specify a space-time connectivity matrix. 
#' 
#' The coefficients \code{beta_ev} are assigned the regularized horseshoe prior (Piironen and Vehtari, 2017), resulting in a relatively sparse model specification. In addition, numerous eigenvectors are automatically dropped because they represent trace amounts of spatial autocorrelation (this is controlled by the \code{threshold} argument). By default, \code{stan_esf} will drop all eigenvectors representing negative spatial autocorrelation patterns. You can change this behavior using the \code{nsa} argument.
#'
#' ### Spatially lagged covariates (SLX)
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
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{data}{a data frame containing the model data}
#' \item{EV}{A matrix of eigenvectors created with \code{w} and \code{geostan::make_EV}}
#' \item{C}{The spatial weights matrix used to construct EV}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re},  the random effects (varying intercepts) formula if provided, and 
#' \code{data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{x_center}{If covariates are centered internally (i.e., `centerx` is not `FALSE`), then `x_centers` is the numeric vector of values on which the covariates were centered.}
#' \item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("esf") and method ("ESF")}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source 
#'
#' Chun, Y., D. A. Griffith, M. Lee and P. Sinha (2016). Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters. *Journal of Geographical Systems*, 18(1), 67-85. \doi{10.1007/s10109-015-0225-3}.
#'
#' Dray, S., P. Legendre & P. R. Peres-Neto (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). *Ecological Modeling*, 196(3-4), 483-493.
#' 
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. *Spatial Statistics*. \doi{10.1016/j.spasta.2020.100450} (open access: \doi{10.31219/osf.io/fah3z}).
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#'
#' Donegan, Connor (2021). Spatial conditional autoregressive models in Stan. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#'
#' Griffith, Daniel A., and P. R. Peres-Neto (2006). Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses. *Ecology* 87(10), 2603-2613.
#' 
#' Griffith, D., and Y. Chun (2014). Spatial autocorrelation and spatial filtering, Handbook of Regional Science. Fischer, MM and Nijkamp, P. eds.
#'
#' Griffith, D., Chun, Y. and Li, B. (2019). *Spatial Regression Analysis Using Eigenvector Spatial Filtering*. Elsevier.
#' 
#' Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In *Electronic Journal of Statistics*, 11(2):5018-5051. \doi{10.1214/17-EJS1337SI}.
#' 
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(sf)
#' library(bayesplot)
#' data(sentencing)
#'
#' # spatial weights matrix with binary coding scheme
#' C <- shape2mat(sentencing, style = "B")
#'
#' # log-expected number of sentences
#' ## expected counts are based on county racial composition and mean sentencing rates
#' log_e <- log(sentencing$expected_sents)
#'
#' # fit spatial Poisson model with ESF + unstructured 'random effects'
#' fit.esf <- stan_esf(sents ~ offset(log_e),
#'                    re = ~ name,
#'                    family = poisson(),
#'                    data = sentencing,
#'                    C = C,
#'                    refresh = 0
#' )
#' 
#' # spatial diagnostics 
#' sp_diag(fit.esf, sentencing)
#'
#' # plot marginal posterior distributions of beta_ev (eigenvector coefficients)
#' plot(fit.esf, pars = "beta_ev")
#'
#' # plot the marginal posterior distributions of the spatial filter (ESF * beta_ev)
#' plot(fit.esf, pars = "beta_ev")
#'
#' # posterior predictive distribution
#' yrep <- posterior_predict(fit.esf, samples = 75)
#' y <- sentencing$sents
#' bayesplot::ppc_dens_overlay(y, yrep) 
#'
#' # map the spatial filter
#' sp.filter <- spatial(fit.esf)$mean
#' st_as_sf(sentencing) %>%
#'  ggplot() +
#'  geom_sf(aes(fill = sp.filter)) +
#'  scale_fill_gradient2()
#'
#' # calculate log-standardized sentencing ratios (log-SSRs)
# (like Standardized Incidence Ratios: observed/exected case counts)
#' f <- fitted(fit.esf)$mean
#' SSR <-  f / sentencing$expected_sents
#' log.SSR <- log( SSR, base = 2 )
#'
#' # map the log-SSRs
#' st_as_sf(sentencing) %>%
#'  ggplot() +
#'  geom_sf(aes(fill = log.SSR)) +
#'  scale_fill_gradient2(
#'    midpoint = 0,
#'    name = NULL,
#'    breaks = seq(-3, 3, by = 0.5)
#'  ) +
#'  labs(title = "Log-Standardized Sentencing Ratios",
#'       subtitle = "log( Fitted/Expected ), base 2"
#'  ) +
#'  theme_void() +
#'  theme(
#'    legend.position = "bottom",
#'    legend.key.height = unit(0.35, "cm"),
#'    legend.key.width = unit(1.5, "cm")
#'  )
#'
#' } 
#'
#' @export
#' @md
#' @importFrom rstan extract_sparse_parts
#' 
stan_esf <- function(formula,
                     slx,
                     re,
                     data,
                     C,
                     EV = make_EV(C, nsa = nsa, threshold = threshold),
                     nsa = FALSE,
                     threshold = 0.25,
                     family = gaussian(),
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
    stopifnot(family$family %in% c("gaussian", "student_t", "poisson", "binomial"))
    stopifnot(!missing(data))
    stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
    stopifnot(all(dim(C) == nrow(data)))
    ## ESF [START] -------------    
    stopifnot(nrow(EV) == nrow(data))
    dev <- ncol(EV)
    ## ESF [STOP] -------------
    a.zero <- as.array(0, dim = 1)
    tmpdf <- as.data.frame(data)
    mod.mat <- model.matrix(formula, tmpdf)
    if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")  
    n <- nrow(mod.mat)
    family_int <- family_2_int(family)
    is_student <- family$family == "student_t"    
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
    ## PRIORS with RHS-ESF [START] -------------
    if (family$family %in% c("poisson", "binomial")) hs_global_scale = 1
    if (family$family %in% c("gaussian", "student_t")) {
        ## default value of global_scale
        p0 <- exp_pars(formula = formula, data = tmpdf, C = C)
        if(p0 >= ncol(EV)) p0 <- ncol(EV) - 0.1
        hs_global_scale <- min(1, p0/(dev-p0) / sqrt(n))
    }
    priors_made <- make_priors(user_priors = prior,
                               y = y,
                               x = x_full,
                               hs_global_scale = hs_global_scale,
                               EV = EV,                               
                               link = family$link,
                               offset = offset)
    standata <- append_priors(standata, priors_made)
    esf_dl <- list(
        dev = dev,
        EV = EV,
        global_scale = priors_made$beta_ev$global_scale,
        slab_df = priors_made$beta_ev$slab_df,
        slab_scale = priors_made$beta_ev$slab_scale
    )
    standata <- c(standata, esf_dl)
    ## PRIORS with RHS-ESF [END] -------------
    ## EMPTY PLACEHOLDERS
    standata <- c(standata, empty_icar_data(n))
    ## ME MODEL STUFF -------------  
    me.list <- prep_me_data(ME, x_no_Wx)
    standata <- c(standata, me.list)
    ## INTEGER OUTCOMES -------------    
    if (family$family == "binomial") {
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
    ## PARAMETERS TO KEEP -------------          
    pars <- c(pars, 'intercept', 'esf', 'beta_ev', 'log_lik', 'fitted')
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
    priors_made_slim <- priors_made[which(names(priors_made) %in% c(pars, "beta_ev"))]
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
    samples <- rstan::sampling(stanmodels$foundation, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
    ## OUTPUT -------------    
    out <- clean_results(samples, pars, is_student, has_re, Wx, x_no_Wx, me.list$x_me_idx)  
    out$data <- ModData
    out$family <- family
    out$formula <- formula
    out$slx <- slx
    out$C <- C
    out$EV <- EV  
    out$re <- re_list
    out$priors <- priors_made_slim
    out$x_center <- attributes(x_full)$`scaled:center`
    out$ME <- list(has_me = me.list$has_me, spatial_me = me.list$spatial_me)
    if (out$ME$has_me) out$ME <- c(out$ME, ME)
    out$spatial <- data.frame(par = "esf", method = "ESF")
    R <- resid(out, summary = FALSE)
    out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C) )    
  return (out)
}


