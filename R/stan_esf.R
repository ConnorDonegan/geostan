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
#' @param re If the model includes a varying intercept specify the grouping variable here using formula syntax, as in \code{~ ID}. The resulting random effects parameter returned is named \code{alpha_re}. The scale parameter for this term is named `alpha_tau`.
#' @param C Spatial connectivity matrix which will be used to calculate eigenvectors, residual spatial autocorrelation as well as any user specified \code{slx} terms or spatial measurement error (ME) models; it will be row-standardized before calculating \code{slx} terms.
#' @param EV A matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial. If provided, still also provide a spatial weights matrix \code{C} for other purposes.  See \link[geostan]{make_EV} and \link[geostan]{shape2mat}.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#'
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list. Errors are assigned a Gaussian probability distribution and the modeled (true) covariate vector is assigned a Student's t model or, if \code{ME$spatial = TRUE}, an auto Gaussian (CAR) model. Elements of the list \code{ME} may include:
#' \describe{
#' 
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one hundred) provide a vector of zeros and ones indicating which columns are bounded. By default the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of any bounded variables.}
#' \item{spatial}{Logical value indicating whether an auto Gaussian (conditional autoregressive (CAR)) model should be used for the covariates (see \link[geostan]{stan_car} for CAR model details). For \code{stan_esf}, this defaults to \code{spatial = TRUE}. If \code{spatial = TRUE} you must provide \code{car_parts} (see below).}
#' \item{car_parts}{A list of data for the CAR model, as returned by [geostan::prep_car_data()] (be sure to return the matrix \code{C}, by using the argument \code{cmat = TRUE}). Only required if \code{spatial=TRUE}.}
#' }
#' 
#' @param nsa Include eigenvectors representing negative spatial autocorrelation? Default \code{nsa = FALSE}. Ignored if \code{EV} is provided.
#' @param threshold Threshold for eigenvector MC value; eigenvectors with values below threshold will be excluded from the candidate set. Default \code{threshold = 0.25}; ignored if \code{EV} is provided. 
#' @param family The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{student_t()} and \code{poisson(link = "log")}, and \code{binomial(link = "logit")}. 
#' @param p0 Number of eigenvector coefficients expected to be far from zero. If this and \code{prior_rhs} are missing, Chun et al.'s (2016) formula will be used to fill this in; see \link[geostan]{exp_pars}. The value of \code{p0} is used to control the prior degree of sparsity in the model.
#'
#' @param prior A named list of parameters for prior distributions. Priors can be set for the intercept (`intercept`), regression coefficients (`beta`), and scale parameter (`sigma`). Users also may set the prior parameters for the degrees of freedom (`nu`) of a Student's t likelihood. For models with so-called random effects, or varying intercept terms, the prior for the scale parameter (`tau`) can be set here. What follows are details on prior models for these parameters:
#' \describe{
#' \item{intercept}{The intercept is assigned a Gaussian prior distribution; provide a numeric vector with location and scale parameters; e.g. \code{prior = list(intercept = c(location = 0, scale = 10))} to assign a Gaussian prior with mean zero and variance 10^2.}
#' 
#' \item{beta}{Regression coefficients are assigned Gaussian prior distributions. Provide a `data.frame` with two columns (location and scale).
#'
#' The order of variables must follow their order of appearance in the model `formula`; e.g., if the formula is `y ~ x1 + x2`, then providing `prior = list(beta = data.frame(location = c(0, 0), scale = c(1, 1))) will assign standard normal prior distributions to both coefficients.
#'
#' Note that if you also use `slx` terms (spatially lagged covariates), then you have to provide priors for them; `slx` terms *precede* other terms in the model formula. For example, providing of `y ~ x1 + x2` with `slx = ~ x1 + x2` is the same as providing `formula = y ~ Wx1 + Wx2 + x1 + x2`, where `Wx1` indicates the spatially lagged covariate.}
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
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. Default \code{chains = 4}.
#' @param iter Number of samples per chain. Default \code{iter = 2000}.
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model. Parameters from the RHS priors include \code{tau} (the global shrinkage parameter) and \code{lambda} (the local shrinkage parameter).
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} except that \code{adapt_delta} is raised to \code{.99} and \code{max_treedepth = 15}.
#' @param silent If \code{TRUE}, suppress printed messages including prior specifications and Stan sampling progress (i.e. \code{refresh=0}). Stan's error and warning messages will still print.
#' @param ... Other arguments passed to \link[rstan]{sampling}. 
#' @details
#'
#' Eigenvector spatial filtering (ESF) is extensivly covered in Griffith et al. (2019). This function implements the methodology introduced in Donegan et al. (2020), drawing on the Piironen and Vehtari's (2017) regularized horseshoe prior.
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
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{data}{a data frame containing the model data}
#' \item{EV}{A matrix of eigenvectors created with \code{w} and \code{geostan::make_EV}}
#' \item{C}{The spatial weights matrix used to construct EV}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ESF component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list containing \code{re},  the random effects (varying intercepts) formula if provided, and 
#' \code{data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
#' \item{spatial}{A data frame with the name of the spatial component parameter ("esf") and method ("ESF")}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#' 
#' @source 
#'
#' Chun, Y., D. A. Griffith, M. Lee and P. Sinha (2016). Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters. *Journal of Geographical Systems*, 18(1), 67-85. \url{10.1007/s10109-015-0225-3}
#'
#' Dray, S., P. Legendre & P. R. Peres-Neto (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). *Ecological Modeling*, 196(3-4), 483-493.
#' 
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. *Spatial Statistics* 38. \url{doi.org/10.1016/j.spasta.2020.100450}
#'
#' Griffith, Daniel A., and P. R. Peres-Neto (2006). Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses. *Ecology* 87(10), 2603-2613.
#' 
#' Griffith, D., and Y. Chun (2014). Spatial autocorrelation and spatial filtering, Handbook of Regional Science. Fischer, MM and Nijkamp, P. eds.
#'
#' Griffith, D., Chun, Y. and Li, B. (2019). *Spatial Regression Analysis Using Eigenvector Spatial Filtering*. Elsevier.
#' 
#' Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In *Electronic Journal of Statistics*, 11(2):5018-5051. \url{https://projecteuclid.org/euclid.ejs/1513306866}
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
stan_esf <- function(formula, slx, re, data, C, EV, ME = NULL,
                     nsa = FALSE, threshold = 0.25,
                     family = gaussian(),
                     p0,
                     prior = NULL,
                     centerx = FALSE, scalex = FALSE,
                     prior_only = FALSE,
                     chains = 4, iter = 2e3, refresh = 500, pars = NULL,
                     control = list(adapt_delta = .99, max_treedepth = 15),
                     silent = FALSE,
                     ...) {
  if (missing(formula)) stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (class(family) != "family" | !family$family %in% c("gaussian", "student_t", "poisson", "binomial")) stop ("Must provide a valid family object: gaussian(), student_t(), or poisson().")
  if (missing(C) | missing(data)) stop ("Must provide data and a spatial connectivity matrix C.")
  if (scalex) centerx <- TRUE
  if (silent) refresh = 0
    ## GLM STUFF -------------  
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")  
  n <- nrow(mod.mat)
  family_int <- family_2_int(family)
  intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0) 
    ## ESF STUFF -------------    
  if (missing(EV)) EV <- make_EV(C, nsa = nsa, threshold = threshold)  
  dev <- ncol(EV)
  if (nrow(EV) != n) stop("nrow(EV) must equal the number of observations.")
    ## GLM STUFF -------------  
  if (intercept_only) {
    if (!missing(slx)) stop("Spatial lag of X (slx) term provided for an intercept only model. Did you intend to include a covariate? If you intend to specify a model in which the only covariate is a spatially-lagged term, you must create this covariate yourself and include it in the main model formula.")      
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
        if (any(rowSums(C) != 1)) {
            message("Creating row-standardized W matrix from C to calculate SLX terms: W = C / rowSums(C)")
        }
        W <- C / rowSums(C)
        Wx <- SLX(f = slx, DF = tmpdf, W = W)
        if (scalex) Wx <- scale(Wx)                     
        dwx <- ncol(Wx)
        dw_nonzero <- sum(W!=0)
        wx_idx <- as.array( which(paste0("w.", dimnames(x)[[2]]) %in% dimnames(Wx)[[2]]), dim = dwx )            
        x <- cbind(Wx, x)
    }
    dbeta_prior <- ncol(x) ## dimensions of beta prior; x includes slx, if any; x.list$x is the processed model matrix without slx terms.
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
  ## ESF STUFF -------------
  if (family$family %in% c("poisson", "binomial")) rhs_scale_global = 1
  if (family$family %in% c("gaussian", "student_t")) {
      # if no prior is provided at all for the rhs:
      if (is.null(prior$rhs)) {
          if (missing(p0)) {
              if (missing(C)) stop("To calculate the prior for the eigenvector coefficients, please provide connectivity matrix C.")
              p0 <- exp_pars(formula = formula, data = tmpdf, C = C)
              if(p0 >= ncol(EV)) p0 <- ncol(EV) - 0.1
              }
       rhs_scale_global <- min(1, p0/(dev-p0) / sqrt(n))
      }
  }
  if (!is.null(prior$rhs)) {
      stopifnot(
          length(prior$rhs) == 3 &
          length(intersect(names(prior$rhs), c("slab_df", "slab_scale", "scale_global"))) == 3
                )
      rhs_scale_global <- as.numeric(prior$rhs["scale_global"])
  }
  ## PARAMETER MODEL STUFF -------------  
  is_student <- family$family == "student_t"
  priors_made <- make_priors(user_priors = prior, y = y, x = x, xcentered = centerx,
                        rhs_scale_global = rhs_scale_global, link = family$link, EV = EV, offset = offset)
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
    alpha_prior = priors_made$intercept,
    beta_prior = t(priors_made$beta), 
    sigma_prior = priors_made$sigma,
    alpha_tau_prior = priors_made$alpha_tau,
    t_nu_prior = priors_made$nu,
    family = family_int,
  ## slx data -------------
    W = W,
    dwx = dwx,
    wx_idx = wx_idx,
  ## esf data -------------
    dev = dev,
    EV = EV,
    slab_df = priors_made$rhs["slab_df"],
    slab_scale = priors_made$rhs["slab_scale"],    
    scale_global = priors_made$rhs["scale_global"],
    prior_only = prior_only
    )
  ## DATA MODEL STUFF -------------  
  me.list <- prep_me_data(ME, x.list$x)
  standata <- c(standata, me.list)
  # -------------  
  # handling multiple possible data types  
  if (family$family == "binomial") {
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'esf', 'beta_ev', 'residual', 'log_lik', 'yrep', 'fitted')
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (family$family %in% c("gaussian", "student_t")) pars <- c(pars, 'sigma')
  if (is_student) pars <- c(pars, "nu")
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  priors_made_slim <- priors_made[which(names(priors_made) %in% c(pars, "rhs"))]
    ## PRINT STUFF -------------    
  if (!silent) print_priors(prior, priors_made_slim)
  ## CALL STAN -------------    
   samples <- rstan::sampling(stanmodels$esf, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  if (missing(C)) C <- NA
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x, me.list$x_me_unbounded_idx, me.list$x_me_bounded_idx)  
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$C <- C
  out$EV <- EV  
  out$re <- re_list
  out$priors <- priors_made_slim
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME
  out$spatial <- data.frame(par = "esf", method = "ESF")
  class(out) <- append("geostan_fit", class(out))
  return (out)
}


