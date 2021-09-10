#' Spatial filtering
#'
#' @export
#'
#' @md
#' 
#' @description Fit a spatial regression model using eigenvector spatial filtering (ESF).
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models are specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{slx = ~ x1 + x2} (the implicit intercept term in this formula is removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept specify the grouping variable here using formula syntax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}. The scale parameter for this term is named `alpha_tau`. That is, \code{alpha_re ~ N(0, alpha_tau)}, \code{alpha_tau ~ Student_t(d.f., location, scale)}.
#' 
#' @param C Spatial connectivity matrix which will be used to calculate eigenvectors, if not provided by the user. Typically, the binary connectivity matrix is best for calculating eigenvectors (i.e., using `C = shape2mat(shape, style = "B")` rather than `style = "W"`). This matrix will also be used to calculate residual spatial autocorrelation (Moran coefficient), any user specified \code{slx} terms, and any spatial measurement error (ME) models; it will be row-standardized before calculating \code{slx} terms. See \code{\link[geostan]{shape2mat}}.
#'
#' @param nsa Include eigenvectors representing negative spatial autocorrelation? Defaults to \code{nsa = FALSE}. This is ignored if \code{EV} is provided.
#' @param threshold Eigenvectors with standardized Moran coefficient values below this `threshold` value will be excluded from the candidate set of eigenvectors, `EV`. This defaults to \code{threshold = 0.25}, and is ignored if \code{EV} is provided. 
#' 
#' @param EV A matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial. With connectivity matrix `C`, you can create `EV` using `EV = make_EV(C)`. If `EV` is provided, still also provide a spatial weights matrix \code{C} for other purposes.  See \code{\link[geostan]{make_EV}} and \code{\link[geostan]{shape2mat}}.
#' 
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'@param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list. Errors are assigned a Gaussian probability distribution and the modeled (true) covariate vector is assigned a Student's t model or, if \code{ME$car_parts} is provided, an auto Gaussian (CAR) model. Elements of the list \code{ME} may include:
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
#' }
#'
#' @param family The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{student_t()} and \code{poisson(link = "log")}, and \code{binomial(link = "logit")}. 
#' @param p0 Number of eigenvector coefficients expected to be far from zero. If this and \code{prior_rhs} are missing, Chun et al.'s (2016) formula will be used to fill this in; see \link[geostan]{exp_pars}. The value of \code{p0} is used to control the prior degree of sparsity in the model.
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
#'
#' \item{rhs}{This is for controlling the regularized horseshoe (RHS) prior, which is assigned to the eigenvector coefficients. Provide a *named* vector of length three. The RHS prior has two parts to consider. First is the 'slab' that is used to regularize large coefficient estimates; this is a zero-mean Student's t model, the user provides the degrees of freedom and scale parameters. The second component is the global shrinkage parameter, controlling the degree of shrinkage; provide a value nearer to zero if you wish to make the model more sparse. To allow the spatial filter to account for a greater amount of spatial autocorrelation (i.e., if you find the residuals contain spatial autocorrelation), increase the global scale parameter up to one. E.g., \code{prior = list(rhs = c(slab_df = 15, slab_scale = 5, scale_global = 0.5))}.}
#' }
#' 
#' @param centerx To center predictors, provide either a logical value (TRUE, FALSE) or numeric-alike vector of length equal to the number of columns of ‘x’, where ‘numeric-alike’ means that ‘as.numeric(.)’ will be applied successfully if ‘is.numeric(.)’ is not true. Passed to \code{\link[base]{scale}}.
#' 
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model. Parameters from the RHS priors include \code{tau} (the global shrinkage parameter) and \code{lambda} (the local shrinkage parameter).
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} except that \code{adapt_delta} is raised to \code{.99} and \code{max_treedepth = 15}.
#' 
#' @param ... Other arguments passed to \link[rstan]{sampling}. 
#' @details
#'
#' Eigenvector spatial filtering (ESF) is a method for spatial regression analysis. ESF is extensivly covered in Griffith et al. (2019). This function implements the methodology introduced in Donegan et al. (2020), which uses Piironen and Vehtari's (2017) regularized horseshoe prior.
#'
#' ESF models take the spectral decomposition of a transformed spatial connectivity matrix, \code{C}. The resulting eigenvectors, `EV`, are mutually orthogonal and uncorrelated map patterns. ESF decomposes spatial autocorrelation into a linear combination of various patterns, typically at different scales (such as local, regional, and global trends). By adding a spatial filter to a regression model, any spatial autocorrelation is shifted from the residuals to the spatial filter. The spatail filter is `EV * beta_ev`, where `beta_ev` is a vector of coefficients.
#'
#' ESF decomposes the data into a global mean, `alpha`, global patterns contributed by covariates, `X * beta`, spatial trends, `EV * beta_ev`, and residual variation. Thus, for `family=gaussian()`,
#' 
#' ```
#'                          Y ~ Gauss(alpha + X * beta + EV * beta_ev, sigma).
#'```
#' An ESF component can be incorporated into the linear predictor of any generalized linear model. For example, a spatial Poisson model for rare disease incidence may be specified as follows:
#' ```
#'                           Y ~ Poisson(exp(offset + Mu))
#'                           Mu = alpha + A + EV * beta_ev
#'                           A ~ Guass(0, tau)
#'                           tau ~ student(20, 0, 2)
#'                           beta_ev ~ horseshoe(.)
#' ```
#' 
#' The \link[geostan]{spatial} method will return `EV * beta`.
#'
#' The model can also be extended to the space-time domain; see \link[geostan]{shape2mat} to specify a space-time connectivity matrix. 
#' 
#' The coefficients \code{beta_ev} are assigned the regularized horseshoe prior (Piironen and Vehtari, 2017), resulting in a relatively sparse model specification. In addition, numerous eigenvectors are automatically dropped because they represent trace amounts of spatial autocorrelation (this is controlled by the \code{threshold} argument). By default, \code{stan_esf} will drop all eigenvectors representing negative spatial autocorrelation patterns. You can change this behavior using the \code{nsa} argument.
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
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. *Spatial Statistics*. \doi{10.1016/j.spasta.2020.100450}.
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
stan_esf <- function(formula,
                     slx,
                     re,
                     data,
                     C,
                     nsa = FALSE,
                     threshold = 0.25,                     
                     EV = make_EV(C, nsa = nsa, threshold = threshold),
                     ME = NULL,
                     family = gaussian(),
                     p0,
                     prior = NULL,
                     centerx = FALSE,
                     prior_only = FALSE,
                     chains = 4, iter = 2e3, refresh = 500, pars = NULL,
                     control = list(adapt_delta = .975, max_treedepth = 13),
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
    ## PRIORS with RHS-ESF [START] -------------
    if (family$family %in% c("poisson", "binomial")) rhs_scale_global = 1
    if (family$family %in% c("gaussian", "student_t")) {
      # if no prior is provided at all for the rhs:
        if (is.null(prior$rhs)) {
            if (missing(p0)) {
                p0 <- exp_pars(formula = formula, data = tmpdf, C = C)
                if(p0 >= ncol(EV)) p0 <- ncol(EV) - 0.1
            }
            rhs_scale_global <- min(1, p0/(dev-p0) / sqrt(n))
        }
    }
    if (!is.null(prior$rhs)) {
        stopifnot(all(c("slab_df", "slab_scale", "scale_global") %in% names(prior$rhs)))
        rhs_scale_global <- as.numeric(prior$rhs["scale_global"])
    }
    is_student <- family$family == "student_t"
    priors_made <- make_priors(user_priors = prior,
                               y = y,
                               x = x_full,
                               rhs_scale_global = rhs_scale_global,
                               EV = EV,                               
                               link = family$link,
                               offset = offset)
    ## PRIORS with RHS-ESF [STOP] -------------
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
    wx_idx = wx_idx,    
  ## esf data -------------
    dev = dev,
    EV = EV,
    slab_df = priors_made$rhs["slab_df"],
    slab_scale = priors_made$rhs["slab_scale"],    
    scale_global = priors_made$rhs["scale_global"]
    )
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
    if (me.list$has_me) pars <- c(pars, "x_true", "mu_x_true", "sigma_x_true")
    if (me.list$spatial_me) pars <- c(pars, "car_rho_x_true")
    priors_made_slim <- priors_made[which(names(priors_made) %in% c(pars, "rhs"))]
    if (me.list$has_me) priors_made_slim <- c(priors_made_slim, list(ME_location = me.list$ME_prior_mean, ME_scale = me.list$ME_prior_scale))    
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


