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
#' alpha_re ~ N(0, alpha_tau)
#' alpha_tau ~ Student_t(d.f., location, scale).
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
#' @param iter Number of samples per chain. .
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples; set \code{refresh=0} to silence this.
#' @param keep_all  If `keep_all = TRUE` then samples for all parameters in the Stan model will be kept; this is necessary if you want to do model comparison with Bayes factors and the `bridgesampling` package.
#' @param slim If `slim = TRUE`, then the Stan model will not collect the most memory-intensive parameters (including n-length vectors of fitted values, log-likelihoods, and ME-modeled covariate values). This will disable many convenience functions that are otherwise available for fitted \code{geostan} models, such as the extraction of residuals, fitted values, and spatial trends, WAIC, and spatial diagnostics, and ME diagnostics; many quantities of interest, such as fitted values and spatial trends, can still be calculated manually using given parameter estimates. The "slim" option is designed for data-intensive routines, such as regression with raster data, Monte Carlo studies, and measurement error models. For more control over which parameters are kept or dropped, use the `drop` argument instead of `slim`.
#' @param drop Provide a vector of character strings to specify the names of any parameters that you do not want MCMC samples for. Dropping parameters in this way can improve sampling speed and reduce memory usage. The following parameter vectors can potentially be dropped from ICAR models:
#' \describe{
#' \item{fitted}{The N-length vector of fitted values}
#' \item{log_lik}{The N-length vector of pointwise log-likelihoods, which is used to calculate WAIC.}
#' \item{alpha_re}{Vector of 'random effects'/varying intercepts.}
#' \item{x_true}{N-length vector of 'latent'/modeled covariate values created for measurement error (ME) models.}
#' \item{phi}{The N-length vector of spatially-autocorrelated parameters (with the ICAR prior).}
#' \item{theta}{The N-length vector of spatially unstructured parameters ('random effects'), for the BYM and BYM2 models.}
#' }
#' If `slim = TRUE`, then `drop` will be ignored.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \code{\link[rstan]{stan}} for details. 
#' @param ... Other arguments passed to \link[rstan]{sampling}.
#' @details
#' 
#'  The intrinsic conditional autoregressive (ICAR) model for spatial data was introduced by Besag et al. (1991). The Stan code for the ICAR component of the model and the BYM2 option is from Morris et al. (2019) with adjustments to enable non-binary weights and disconnected graph structures (see Freni-Sterrantino (2018) and Donegan (2021)).
#'
#' The exact specification depends on the `type` argument. 
#'
#' ### ICAR
#'
#' For Poisson models for count data, y, the basic model specification (`type = "icar"`) is:
#' 
#' \deqn{y ~ Poisson(e^{O + \mu + \phi}) }
#' \deqn{\phi \sim ICAR(\tau_s) }
#' \deqn{\tau_s \sim Gauss(0, 1)}
#' 
#'  where \eqn{\mu} contains an intercept and potentially covariates. The spatial trend \eqn{phi} has a mean of zero and a single scale parameter \eqn{\tau_s} (which user's will see printed as the parameter named `spatial_scale`).
#' 
#' The ICAR prior model is a CAR model that has a spatial autocorrelation parameter \eqn{\rho} equal to 1 (see \link[geostan]{stan_car}). Thus the ICAR prior places high probability on a very smooth spatially (or temporally) varying mean. This is rarely sufficient to model the amount of variation present in social and health data. For this reason, the BYM model is typically employed.
#'
#' ### BYM
#'
#' Often, an observational-level random effect term, `theta`, is added to capture (heterogeneous or unstructured) deviations from \eqn{\mu + \phi}. The combined term is referred to as a convolution term:
#' 
#' \eqn{
#'  convolution = \phi + \theta.
#' }
#' 
#' This is known as the BYM model (Besag et al. 1991), and can be specified using `type = "bym"`:
#'
#' \eqn{y \sim Poisson(e^{O + \mu + \phi + \theta}) }
#' \deqn{ \phi \sim ICAR(\tau_s) }
#' \deqn{ \theta \sim Gaussian(0, \tau_{ns}) }
#' \deqn{ \tau_s \sim Gaussian(0, 1) }
#' \deqn{ \tau_{ns} \sim Gaussian(0, 1) }
#'
#' The model is named after Besag, York, and Mollié (1991).
#' 
#' ### BYM2
#' 
#' Riebler et al. (2016) introduce a variation on the BYM model (`type = "bym2"`). This specification combines \eqn{\phi} and \eqn{\theta} using a mixing parameter \eqn{\rho} that controls the proportion of the variation that is attributable to the spatially autocorrelated term \eqn{\phi} rather than the spatially unstructured term \eqn{\theta}. The terms share a single scale parameter \eqn{\tau}:
#' 
#' \deqn{convolution = [sqrt(\rho * S) * \tilde{\phi} + sqrt(1 - \rho)  \tilde{\theta}] * \tau
#' }
#' \deqn{ \tilde{\phi} \sim Gaussian(0, 1) }
#' \deqn{ \tilde{\theta} \sim Gaussian(0, 1) }
#' \deqn{ \tau \sim Gaussian(0, 1) }
#' 
#' The terms \eqn{\tilde{\phi}}, \eqn{\tilde{\theta}} are standard normal deviates, \eqn{\rho} is restricted to values between zero and one, and \eqn{S} is the 'scale_factor' (a constant term provided by the user). By default, the 'scale_factor' is equal to one, so that it does nothing. Riebler et al. (2016) argue that the interpretation or meaning of the scale of the ICAR model depends on the graph structure of the connectivity matrix \eqn{C}. This implies that the same prior distribution assigned to \eqn{\tau_s} will differ in its implications if \eqn{C} is changed; in other words, the priors are not transportable across models, and models that use the same nominal prior actually have different priors assigned to \eqn{\tau_s}.
#'
#' Borrowing `R` code from Morris (2017) and following Freni-Sterrantino et al. (2018), the following `R` code can be used to create the 'scale_factor' \eqn{S} for the BYM2 model (note, this requires the INLA R package), given a spatial adjacency matrix, \eqn{C}:
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
#' This code adjusts for 'islands' or areas with zero neighbors, and it also handles disconnected graph structures (see Donegan and Morris 2021). Following Freni-Sterrantino (2018), disconnected components of the graph structure are given their own intercept term; however, this value is added to \eqn{\phi} automatically inside the Stan model. Therefore, the user never needs to make any adjustments for this term. (To avoid complications from using a disconnected graph structure, you can apply a proper CAR model instead of the ICAR: \code{\link[geostan]{stan_car}}).
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
#' ## Additional functionality
#'
#' The CAR models can also incorporate spatially-lagged covariates, measurement/sampling error in covariates (particularly when using small area survey estimates as covariates), and censored outcomes (such as arise when a disease surveillance system suppresses data for privacy reasons). For details on these options, please see the Details section in the documentation for \link[geostan]{stan_glm}.
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors and the weight attached to each pair (i.e., their corresponding element in the connectivity matrix  C}
#' \item{C}{Spatial connectivity matrix}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ICAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list with two name elements, \code{formula} and \code{Data}, containing the formula \code{re} and a data frame with columns \code{id} (the grouping variable) and \code{idx} (the index values assigned to each group).}
#' \item{priors}{Prior specifications.}
#' \item{x_center}{If covariates are centered internally (`centerx = TRUE`), then `x_center` is a numeric vector of the values on which covariates were centered.}
#' \item{spatial}{A data frame with the name of the spatial parameter (\code{"phi"} if \code{type = "icar"} else \code{"convolution"}) and method (\code{toupper(type)}).}
#' }
#' 
#' @author Connor Donegan, \email{connor.donegan@gmail.com}
#'
#' @seealso \link[geostan]{shape2mat}, \link[geostan]{stan_car}, \link[geostan]{stan_esf}, \link[geostan]{stan_glm}, \link[geostan]{prep_icar_data}
#' 
#' @source
#'
#' Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. *Journal of the Royal Statistical Society: Series B (Methodological)*, 36(2), 192-225.
#'
#' Besag, J., York, J., and Mollié, A. (1991). Bayesian image restoration, with two applications in spatial statistics. *Annals of the Institute of Statistical Mathematics*, 43(1), 1-20.
#'
#' Donegan, Connor and Morris, Mitzi (2021). Flexible functions for ICAR, BYM, and BYM2 models in Stan. Code repository. <https://github.com/ConnorDonegan/Stan-IAR>
#'
#' Donegan, Connor (2021b). Building spatial conditional autoregressive (CAR) models in the Stan programming language. *OSF Preprints*. \doi{10.31219/osf.io/3ey65}.
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#' 
#' Freni-Sterrantino, Anna, Massimo Ventrucci, and Håvard Rue (2018). A Note on Intrinsic Conditional Autoregressive Models for Disconnected Graphs. *Spatial and Spatio-Temporal Epidemiology*, 26: 25–34.
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. *Spatial and spatio-temporal epidemiology*, 31, 100301.
#'
#' Riebler, A., Sorbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. *Statistical Methods in Medical Research*, 25(4), 1145-1165.
#'
#' @examples
#' \donttest{
#' data(sentencing)
#' C <- shape2mat(sentencing, "B")
#' log_e <- log(sentencing$expected_sents)
#' fit.bym <- stan_icar(sents ~ offset(log_e),
#'                      family = poisson(),
#'                      data = sentencing,
#'                      type = "bym",
#'                      C = C,
#'                      chains = 2, iter = 800) # for speed only
#'
#' # spatial diagnostics
#' sp_diag(fit.bym, sentencing)
#'                                        
#' # check effective sample size and convergence
#' library(rstan)
#' rstan::stan_ess(fit.bym$stanfit)
#' rstan::stan_rhat(fit.bym$stanfit)
#' 
#' # calculate log-standardized incidence ratios 
#' # (observed/exected case counts)
#' library(ggplot2)
#' library(sf)
#' 
#' f <- fitted(fit.bym, rates = FALSE)$mean
#' SSR <- f / sentencing$expected_sents
#' log.SSR <- log( SSR, base = 2)
#'
#' ggplot(sentencing) +
#'   geom_sf(aes(fill = log.SSR)) +
#'   scale_fill_gradient2(
#'    low = "navy",
#'    high = "darkred"
#'   ) +
#'   labs(title = "Log-standardized sentencing ratios",
#'        subtitle = "log( Fitted/Expected), base 2") +
#'   theme_void() +
#'   theme(
#'    legend.position = "bottom",
#'    legend.key.height = unit(0.35, "cm"),
#'    legend.key.width = unit(1.5, "cm")
#'   )
#' }
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
                      censor_point,
                      prior_only = FALSE,
                      chains = 4, iter = 2e3, refresh = 500,
                      keep_all = FALSE,
                      slim = FALSE,
                      drop = NULL,
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
        dw_nonzero = length(W.list$w),
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
    ## ICAR DATA [START] -------------
    iar.list <- prep_icar_data(C, scale_factor = scale_factor)
    if (!inherits(scale_factor, "NULL")) stopifnot(length(scale_factor) == iar.list$k)
    standata <- c(standata, iar.list)
    standata$type <- match(type, c("icar", "bym", "bym2"))
    ## ICAR DATA [STOP] -------------    
    ## EMPTY PLACEHOLDERS
    standata <- c(standata, empty_esf_data(n), empty_car_data(), empty_sar_data(n))
    ## ME MODEL -------------  
    me.list <- make_me_data(ME, xraw)
    standata <- c(standata, me.list)  
    ## PARAMETERS TO KEEP with ICAR [START] -------------        
    pars <- c(pars, 'intercept', 'spatial_scale', 'fitted', 'phi', 'log_lik')
    if (type == "bym2") pars <- c(pars, 'theta', 'rho')
    if (type == "bym") pars <- c(pars, 'theta', 'theta_scale')
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
    if (slim == TRUE) drop <- c('fitted', 'log_lik', 'alpha_re', 'x_true', 'phi', 'theta')
    pars <- drop_params(pars = pars, drop_list = drop)
    priors_made_slim <- priors_made[which(names(priors_made) %in% pars)]
    if (me.list$has_me) priors_made_slim$ME_model <- ME$prior        
    print_priors(prior, priors_made_slim)    
    ## PARAMETERS TO KEEP with ICAR [STOP] -------------
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
    out$C <- as(C, "dMatrix")    
    out$re <- re_list
    out$priors <- priors_made_slim
    out$x_center <- get_x_center(standata, samples)
    out$ME <- list(has_me = me.list$has_me, spatial_me = me.list$spatial_me)
    if (out$ME$has_me) out$ME <- c(out$ME, ME)
    ## ICAR OUTPUT [START] --------
    out$edges <- edges(C)       
    out$spatial <- data.frame(par = "phi", method = toupper(type))
    ## ICAR OUTPUT [STOP] --------
    if (!missing(C) && any(pars == 'fitted')) {
        R <- resid(out, summary = FALSE)
        out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C, warn = FALSE, na.rm = TRUE) )
    }
    out$N <- length( y_index_list$y_obs_idx )    
    return (out)
}

