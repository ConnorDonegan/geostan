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
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally); you must also provide \code{C} when including \code{slx}.
#'  Specified covariates will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term (or "spatially unstructured random effect") specify the grouping variable here using formula syntax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}, and the associated scale parameter is \code{alpha_tua}. That is, \code{alpha_re ~ N(0, alpha_tau)}, \code{alpha_tau ~ Student_t(d.f., location, scale)}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list. Errors are assigned a Gaussian probability distribution and the modeled (true) covariate vector is assigned a Student's t model or, if \code{ME$car_parts} is provided, an auto Gaussian (CAR) model. Elements of the list \code{ME} may include:
#' \describe{
#' 
#' \item{se}{A dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' 
#' \item{bounds}{An optional numeric vector of length two providing the upper and lower bounds, respectively, of the variables. If not provided, they will be set to `c(-Inf, Inf)` (i.e., unbounded). Common usages include keeping percentages between zero and one hundred or proportions between zero and one.}
#' 
#'  \item{prior}{Provide parameter values for the prior distributions of the measurement error model(s). The spatial conditional autoregressive (CAR) and non-spatial student's t models both contain a location parameter (the mean) and a scale parameter, each of which require prior distributions. If none are provided, default priors will be assigned and printed to the console. 
#'
#' To set these priors to custom values, provide a named list with items `location` and `scale` (you must provide both). The prior for the location parameter is Gaussian, the default being `Gauss(0, 100)`. To change these values, provide a `data.frame` with columns named `location` and `scale`. Provide prior specifications for each covariate in `se`; list the priors in the same order as the columns of `se`.
#'
#' The prior for the `scale` parameters is Student's t, and the default is `Student_t(10, 0, 40)`. The degrees of freedom (10) and mean (zero) are fixed, but you can alter the scale by providing a vector of values in the same order as the columns of `se`.
#'
#' For example, `ME$prior <- list(location = data.frame(location = c(0, 0), scale = c(100, 100))); ME$prior$scale <- c(40, 40)`.
#'
#' The CAR model also has a spatial autocorrelation parameter, `car_rho`, which is assigned a uniform prior distribution. You can set the boudaries of the prior with `ME$prior$car_rho <- c(lower_bound, upper_bound)`. 
#' }
#' 
#' \item{car_parts}{A list of data for the CAR model, as returned by \link[geostan]{prep_car_data}. If not provided, a non-spatial Student's t model will be used instead of the CAR model.}
#'
#' }
#' @param C Optional spatial connectivity matrix which will be used to calculate residual spatial autocorrelation as well as any user specified \code{slx} or spatial measurement error (\code{ME}) terms; it will automatically be row-standardized before calculating \code{slx} terms.  See \code{\link[geostan]{shape2mat}}.
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{poisson(link = "log")}, \code{binomial(link = "logit")}, \code{student_t()}, and the default \code{gaussian()}.
#' 
#' @param prior A named list of parameters for prior distributions. User-defined priors can be assigned to the following parameters:
#' \describe{
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution; provide a numeric vector with location and scale parameters; e.g. \code{prior = list(intercept = c(location = 0, scale = 10))} to assign a Gaussian prior with mean zero and variance 10^2.}
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Provide a `data.frame` with two columns (location and scale).
#'
#' The order of variables must follow their order of appearance in the model `formula`; e.g., if the formula is `y ~ x1 + x2`, then providing `prior = list(beta = data.frame(location = c(0, 0), scale = c(1, 1))) will assign standard normal prior distributions to both coefficients.
#'
#' Note that if you also use `slx` terms (spatially lagged covariates), then you have to provide priors for them; `slx` terms *precede* other terms in the model formula. For example, providing of `y ~ x1 + x2` with `slx = ~ x1 + x2` is the same as providing `formula = y ~ I(W \%*\% x1) + I(W \%*\% x2) + x1 + x2`, where `W` is a row-standardized spatial weights matrix (see \code{\link[geostan]{shape2mat}}).}
#'
#' \item{sigma}{The scale parameter for Gaussian and Student's t likelihood models is assigned a half-Student's t prior distribution. Thus, provide values for the degrees of freedom, location, and scale parameters; e.g., `prior = list(sigma = c(df = 10, location = 0, scale = 5))` for a prior centered on zero with scale of 5 and 10 degrees of freedom. The half-Student's t prior for `sigma` is constrained to be positive.}
#'
#' \item{nu}{The degrees of freedom for the Student's t likelihood model is assigned a gamma prior distribution. The default prior is `prior = list(nu = c(alpha = 3, beta = 0.2))`.}
#'
#' \item{tau}{The scale parameter for random effects, or varying intercepts, terms. This scale parameter, `tau`, is assigned a half-Student's t prior. To set this, use, e.g., `prior = list(tau = c(df = 20, location = 0, scale = 20))`.}
#' }
#' 
#' @param centerx To center predictors, provide either a logical value (TRUE, FALSE) or numeric-alike vector of length equal to the number of columns of ‘x’, where ‘numeric-alike’ means that ‘as.numeric(.)’ will be applied successfully if ‘is.numeric(.)’ is not true. Passed to \code{\link[base]{scale}}.
#' 
#' @param prior_only Draw samples from the prior distributions of parameters only. 
#' @param chains Number of MCMC chains to estimate. 
#' @param iter Number of samples per chain. 
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} except that \code{adapt_delta} is raised to \code{0.95} and \code{max_treedepth = 15}.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' 
#' @details
#'
#' Fit a generalized linear model using the R formula interface. Default prior distributions are designed to be weakly informative relative to the data. Much of the functionality intended for spatial models, such as the ability to add spatially lagged covariates and observational error models, are also available in \code{stan_glm}. All of \code{geostan}'s spatial models build on top of the same Stan code used in \code{stan_glm}.
#'
#' Poisson models are often used to calculate incidence rates (mortality rates, or disease incidence rates) for administrative areas, like counties or census tracts; the user provided offset should be, in that case, the natural logarithm of the denominator, e.g., log-population at risk. If `Y` are counts of cases, and `P` are populations at risk, then the crude rates are `Y/P`, and the user should provide an offset value of `log(P)`. For such a case, disease incidence across the collection of areas could be modeled as:
#' ```
#'                           Y ~ Poisson(exp(offset + Mu))
#'                           Mu = alpha + A
#'                           A ~ Guass(0, tau)
#'                           tau ~ student(20, 0, 2),
#' ```
#' where `alpha` is the mean log-risk (incidence rate) and `A` is a vector of (so-called) random effects, which enable partial pooling of information across observations. See the example section of this documenet for a demonstration. 
#'
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re}, the random effects (varying intercepts) formula if provided, and 
#'  \code{Data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{x_center}{If covariates are centered internally (i.e., `centerx` is not `FALSE`), then `x_centers` is the numeric vector of values on which the covariates were centered.}
#' \item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
#' \item{spatial}{NA, slot is maintained for use in \code{geostan_fit} methods.}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(sf)
#' data(sentencing)
#'
#' sentencing$log_e <- log(sentencing$expected_sents)
#' fit.pois <- stan_glm(sents ~ offset(log_e),
#'                      re = ~ name,
#'                      family = poisson(),
#'                      data = sentencing
#'  )
#'
#' # MCMC diagnostics plot: Rhat values should all by very near 1
#' rstan::stan_rhat(fit.pois$stanfit)
#'  # see effective sample size for all parameters and generated quantities
#'  # (including residuals, predicted values, etc.)
#' rstan::stan_ess(fit.pois$stanfit)
#' # or for a particular parameter
#' rstan::stan_ess(fit.pois$stanfit, "alpha_re")
#'
#' # Spatial autocorrelation/residual diagnostics
#' sp_diag(fit.pois, sentencing)
#' 
#' ## Posterior predictive check                                        
#' library(bayesplot)
#' yrep <- posterior_predict(fit.pois, samples = 75)
#' y <- sentencing$sents
#' ppc_dens_overlay(y, yrep)
#' }
stan_glm <- function(formula,
                     slx,
                     re,
                     data,
                     C,
                     ME = NULL,                     
                     family = gaussian(),
                     prior = NULL,
                     centerx = FALSE, 
                     prior_only = FALSE,
                     chains = 4,
                     iter = 2e3, 
                     refresh = 1e3,
                     pars = NULL,
                     control = list(adapt_delta = 0.95, max_treedepth = 12),
                     ...) {
    stopifnot(inherits(formula, "formula"))
    stopifnot(inherits(family, "family"))
    stopifnot(family$family %in% c("gaussian", "student_t", "poisson", "binomial"))
    stopifnot(!missing(data))
    if (!missing(C)) {
        stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
        stopifnot(all(dim(C) == nrow(data)))
    }
    if (!missing(ME)) if (inherits(ME$car_parts$C, "Matrix") | inherits(ME$car_parts$C, "matrix")) C <- ME$car_parts$C
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
        dbeta_prior <- 0
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
      dbeta_prior <- ncol(x_full) 
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
    ## PRIORS -------------  
    is_student <- family$family == "student_t"
    priors_made <- make_priors(user_priors = prior,
                               y = y,
                               x = x_full,
                               link = family$link,
                               offset = offset)
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
    alpha_prior = priors_made$intercept,
    beta_prior = t(priors_made$beta), 
    sigma_prior = priors_made$sigma,
    alpha_tau_prior = priors_made$alpha_tau,
    t_nu_prior = priors_made$nu,
    family = family_int,
    prior_only = prior_only,    
  ## slx data -------------    
    W_w = as.array(W.list$w),
    W_v = as.array(W.list$v),
    W_u = as.array(W.list$u),
    dw_nonzero = length(W.list$w),
    dwx = dwx,
    wx_idx = wx_idx
    )
    ## EMPTY PLACEHOLDERS
    standata <- c(standata, empty_icar_data(n), empty_esf_data(n))
    ## ME MODEL -------------  
    me.list <- prep_me_data(ME, x_no_Wx)
    standata <- c(standata, me.list)  
    ## INTEGER OUTCOMES -------------    
    if (family$family == "binomial") {
        standata$y <- standata$y_int <- y[,1]
        standata$trials <- y[,1] + y[,2]
    }
    ## PARAMETERS TO KEEP -------------               
    pars <- c(pars, 'intercept', 'log_lik', 'fitted')
    if (!intercept_only) pars <- c(pars, 'beta')
    if (dwx) pars <- c(pars, 'gamma')
    if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
    if (is_student) pars <- c(pars, "nu")
    if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
    if (me.list$has_me) pars <- c(pars, "x_true", "mu_x_true", "sigma_x_true")
    if (me.list$spatial_me) pars <- c(pars, "car_rho_x_true")
    priors_made_slim <- priors_made[which(names(priors_made) %in% pars)]
    if (me.list$has_me) priors_made_slim <- c(priors_made_slim, list(ME_location = me.list$ME_prior_mean, ME_scale = me.list$ME_prior_scale))
    if (me.list$has_me && me.list$spatial_me) priors_made_slim <- c(priors_made_slim, list(ME_car_rho = me.list$ME_prior_car_rho))
    print_priors(prior, priors_made_slim)
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
    if (has_re) {
        out$spatial <- data.frame(par = "alpha_re", method = "Exchangeable")
    } else {
        out$spatial <- data.frame(par = "none", method = "none")
    }
    if (!missing(C) && (inherits(C, "matrix") | inherits(C, "Matrix"))) {
        R <- resid(out, summary = FALSE)
        out$diagnostic["Residual_MC"] <- mean( apply(R, 1, mc, w = C) )
    }        
    return(out)
}

