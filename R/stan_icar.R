#' Intrinsic autoregressive models
#'
#' @export
#'
#' @md
#' 
#' @description Assign the intrinsic conditional auto-regressive (ICAR) prior model to parameters. Options include the BYM model, the BYM2 model, and a solo ICAR term. 
#' 
#' @param formula A model formula, following the R \link[stats]{formula} syntax. Binomial models can be specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.
#' @param slx Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally).
#'  These will be pre-multiplied by a row-standardized spatial weights matrix and then added (prepended) to the design matrix.
#'  If and when setting priors for \code{beta} manually, remember to include priors for any SLX terms as well.
#' @param re If the model includes a varying intercept term \code{alpha_re} specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re ~ N(0, alpha_tau)}, \code{alpha_tau ~ Student_t(d.f., location, scale)}. Before using this, read the \code{Details} section and the \code{type} argument.
#' @param data A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.
#' @param type Defaults to "icar" (partial pooling of neighboring observations through parameter \code{phi}); specify "bym" to add a second parameter vector \code{theta} to perform partial pooling across all observations; specify "bym2" for the innovation introduced by Riebler et al. (2016). See \code{Details} for more information.
#' @param scale_factor For the BYM2 model, optional. If missing, this will be set to a vector of ones. 
#' @param ME To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a named list. Errors are assigned a Gaussian probability distribution and the modeled (true) covariate vector is assigned a Student's t model or, if \code{ME$spatial = TRUE}, an auto Gaussian (CAR) model. Elements of the list \code{ME} may include:
#' \describe{
#' 
#' \item{se}{a dataframe with standard errors for each observation; columns will be matched to the variables by column names. The names should match those from the output of \code{model.matrix(formula, data)}.}
#' \item{bounded}{If any variables in \code{se} are bounded within some range (e.g. percentages ranging from zero to one hundred) provide a vector of zeros and ones indicating which columns are bounded. By default the lower bound will be 0 and the upper bound 100, for percentages.}
#' \item{bounds}{A numeric vector of length two providing the upper and lower bounds, respectively, of any bounded variables.}
#' \item{spatial}{Logical value indicating whether an auto Gaussian (i.e., conditional autoregressive (CAR)) model should be used for the covariates. For \code{stan_icar}, this defaults to \code{spatial = TRUE}; if \code{spatial = TRUE} you must provide \code{car_parts} (see below).}
#' \item{car_parts}{A list of data for the CAR model, as returned by \link[geostan]{prep_car_data} (be sure to return the matrix \code{C}, by using the argument \code{cmat = TRUE}). Only required if \code{spatial=TRUE}.}
#' }
#' @param C Spatial connectivity matrix which will be used to construct an edge list for the ICAR model, and to calculate residual spatial autocorrelation as well as any user specified \code{slx} terms or spatial measurement error (ME) models. It will automatically be row-standardized before calculating \code{slx} terms. \code{C} must be a binary symmetric \code{n x n} matrix.
#' 
#' @param family The likelihood function for the outcome variable. Current options are \code{binomial(link = "logit")} and \code{poisson(link = "log")}. 
#' @param prior A \code{data.frame} or \code{matrix} with location and scale parameters for Gaussian prior distributions on the model coefficients. Provide two columns---location and scale---and a row for each variable in their order of appearance in the model formula. Default priors are weakly informative relative to the scale of the data.
#' @param prior_intercept A vector with location and scale parameters for a Gaussian prior distribution on the intercept; e.g. \code{prior_intercept = c(0, 10)}. 
#' @param prior_tau Set hyperparameters for the scale parameter of exchangeable random effects/varying intercepts. The random effects are given a normal prior with scale parameter \code{alpha_tau}. The latter is given a half-Student's t prior with default of 20 degrees of freedom, centered on zero and scaled to the data to be weakly informative. To adjust it use, e.g., \code{prior_tau = c(df = 15, location = 0, scale = 5)}.
#' @param centerx Should the covariates be centered prior to fitting the model? Defaults to \code{FALSE}.
#' @param scalex Should the covariates be centered and scaled (divided by their standard deviation)? Defaults to \code{FALSE}.
#' @param prior_only Draw samples from the prior distributions of parameters only.
#' @param chains Number of MCMC chains to estimate. 
#' @param iter Number of samples per chain. .
#' @param refresh Stan will print the progress of the sampler every \code{refresh} number of samples; set \code{refresh=0} to silence this.
#' @param pars Optional; specify any additional parameters you'd like stored from the Stan model.
#' @param control A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details. The defaults are the same \code{rstan::stan} except that \code{adapt_delta} is raised to \code{.9} and \code{max_treedepth = 15}.
#' @param silent If \code{TRUE}, suppress printed messages including prior specifications and Stan sampling progress (i.e. \code{refresh=0}). Stan's error and warning messages will still print.
#' @param ... Other arguments passed to \link[rstan]{sampling}. For multi-core processing, you can use \code{cores = parallel::detectCores()}, or run \code{options(mc.cores = parallel::detectCores())} first.
#' @details
#' 
#'  The Stan code for the ICAR component of the model and the BYM2 option is from Morris et al. (2019) with adjustments to enable non-binary weights and disconnected graph structures (see Freni-Sterrantino (2018) and Donegan (2021)).
#'
#' The exact specification depends on the `type` argument.
#'
#' ### `type = 'icar'`
#'
#' For Poisson models for count data, y, the basic model specification (`type = "icar"`) is:
#' ```
#'              y ~ Poisson(exp(offset + mu + phi))
#'              phi ~ ICAR(spatial_scale)
#'              spatial_scale ~ Gaussian(0, 1)
#' ```
#'  where `mu` contains an intercept and potentially covariates. The spatial trend, `phi`, has a mean of zero and a single scale parameter, `spatial_scale`.
#' 
#' The ICAR prior model is a CAR model that has a spatial autocorrelation parameter \code{car_alpha} equal to 1 (see \link[geostan]{stan_car}). Thus the ICAR prior places high probability on a smooth spatially (or temporally) varying mean.
#'
#' ### `type = 'bym'`
#'
#' Often, an observational-level random effect term, `theta`, is added to capture (heterogeneous or unstructured) deviations from `mu + phi`. The combined term is referred to as a convolution term:
#' ```
#'              convolution = phi + theta.
#' ```
#' This is known as the BYM model (Besag et al. 1991), and can be specified using `type = "bym"`:
#' ```
#'              y ~ Poisson(exp(offset + mu + phi + theta))
#'              phi ~ ICAR(spatial_scale)
#'              theta ~ Gaussian(0, theta_scale)
#'              spatial_scale ~ Gaussian(0, 1)
#'              theta_scale ~ Gaussian(0, 1)
#' ```
#'
#' ### `type = 'bym2'`
#' 
#' Riebler et al. (2016) introduce a variation on the BYM model (`type = "bym2"`). This specification combines `phi` and `theta` using a mixing parameter, `rho`, that controls the proportion of the variation that is attributable to the spatially autocorrelated term, `phi`, rather than the spatially unstructured term, `theta`. The terms share a single scale parameter:
#' ```
#'              convolution = [sqrt(rho/scale_factor) * phi_tilde + sqrt(1 - rho) * theta_tilde] * spatial_scale.
#'              phi_tilde ~ Gaussian(0, 1)
#'              theta_tilde ~ Gaussian(0, 1)
#' ```
#' The two `_tilde` terms are equivalent to standard normal deviates, `rho` is restricted to values between zero and one, and `scale_factor` is a constant term provided by the user. By default `scale_factor` is equal to one, so that it does nothing. Riebler et al. (2016) argue that the interpretation or meaning of the scale of the ICAR model depends on the graph structure, `C`. This implies that the same prior distribution assigned to the `spatial_scale` will differ in its implications if `C` is changed; in other words, the priors are not transportable across models, and models that use the same nominal prior actually have different priors assigned to `spatial_scale`.
#'
#' Borrowing `R` code from Morris (2017) and following Freni-Sterrantino et al. (2018), the following `R` code can be used to create the `scale_factor` for the BYM2 model (note, this requires the INLA R package):
#' ```
#'               ## create a list of data for stan_icar
#'               icar.data <- geostan::prep_icar_data(C)
#'               ## calculate scale_factor for each of k connected group of nodes, using the scale_c function (Morris et al. 2019)
#'               k <- icar.data$k
#'               scale_factor <- vector(mode = "numeric", length = k)
#'               for (j in 1:k) {
#'                 g.idx <- which(icar.data$comp_id == j) 
#'                 if (length(g.idx) == 1) {
#'                   scale_factor[j] <- 1
#'                   next
#'                  }    
#'               Cg <- C[g.idx, g.idx] 
#'               scale_factor[j] <- scale_c(Cg) 
#'}
#'                ## update the data list for stan_icar: exactly like this
#'               icar.data$inv_sqrt_scale_factor <- 1 / sqrt( scale_factor )
#' ```
#' This code adjusts for 'islands' or areas with zero neighbors, and it also handles disconnected graph structures (see Donegan 2021). Following Freni-Sterrantino (2018), disconnected components of the graph structure are given their own intercept term; however, this value is added to `phi` automatically inside the Stan model. Therefore, the use never needs to make any adjustments for this term.
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
#' 
#' #'   Morris, Mitzi (2017). Spatial Models in Stan: Intrinsic Auto-Regressive Models for Areal Data. <https://mc-stan.org/users/documentation/case-studies/icar_stan.html>
#' #'
#' scale_c <- function(C) {
#'  #' compute geometric mean of a vector
#'  geometric_mean <- function(x) exp(mean(log(x))) 
#'  N = dim(C)[1]
#'  # Create ICAR precision matrix  (diag - C): this is singular
#'  # function Diagonal creates a square matrix with given diagonal
#'  Q =  Diagonal(N, rowSums(C)) - C
#'  # Add a small jitter to the diagonal for numerical stability (optional but recommended)
#'  Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
#'  # Function inla.qinv provides efficient way to calculate the elements of the
#'  # the inverse corresponding to the non-zero elements of Q
#'  Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))
#'  # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
#'  scaling_factor <- geometric_mean(Matrix::diag(Q_inv)) 
#'  return(scaling_factor) 
#'}
#' ```
#' 
#' @return An object of class class \code{geostan_fit} (a list) containing: 
#' \describe{
#' \item{summary}{Summaries of the main parameters of interest; a data frame}
#' \item{diagnostic}{Widely Applicable Information Criteria (WAIC) with crude measure of effective number of parameters (\code{eff_pars}) and 
#'  mean log pointwise predictive density (\code{lpd}), and residual spatial autocorrelation (Moran coefficient of the residuals). Residuals are relative to the mean posterior fitted values.}
#' \item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
#' \item{data}{a data frame containing the model data}
#' \item{edges}{The edge list representing all unique sets of neighbors and the weight attached to each pair (i.e., their corresponding element in the connectivity matrix  C}
#' \item{family}{the user-provided or default \code{family} argument used to fit the model}
#' \item{formula}{The model formula provided by the user (not including ICAR component)}
#' \item{slx}{The \code{slx} formula}
#' \item{re}{A list with two name elements, \code{formula} and \code{Data}, containing the formula \code{re} and a data frame with columns \code{id} (the grouping variable) and \code{idx} (the index values assigned to each group).}
#' \item{priors}{Prior specifications.}
#' \item{scale_params}{A list with the center and scale parameters returned from the call to \code{base::scale} on the model matrix. If \code{centerx = FALSE} and \code{scalex = FALSE} then it is an empty list.}
#' \item{spatial}{A data frame with the name of the spatial parameter (\code{"phi"} if \code{type = "icar"} else \code{"convolution"}) and method (\code{toupper(type)}).}
#' }
#' 
#' @author Connor Donegan, \email{Connor.Donegan@UTDallas.edu}
#'
#' @seealso \link[geostan]{prep_icar_data}, \link[geostan]{shape2mat}, \link[geostan]{stan_car}, \link[geostan]{stan_esf}, \link[geostan]{stan_glm}
#' 
#' @source
#'
#' Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. Journal of the Royal Statistical Society: Series B (Methodological), 36(2), 192-225.
#'
#' Besag, J., York, J., & Mollié, A. (1991). Bayesian image restoration, with two applications in spatial statistics. Annals of the institute of statistical mathematics, 43(1), 1-20.
#'
#' Donegan, Connor. 2021. Flexible functions for ICAR, BYM, and BYM2 models in Stan. Code repository. <https://github.com/ConnorDonegan/Stan-IAR>
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
#'                      C = C,
#'                      refresh = 0
#'  )
#'
#' # check effective sample size and convergence
#' rstan::stan_ess(fit.bym$stanfit)
#' rstan::stan_rhat(fit.bym$stanfit)
#'
#' # see some spatial diagnostics
#' sp_diag(fit.bym, sentencing)
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
stan_icar <- function(formula, slx, re, data,
                      type = c("icar", "bym", "bym2"),
                      scale_factor = NULL,
                      ME = NULL, C, 
                      family = poisson(),
                      prior = NULL, prior_intercept = NULL, prior_tau = NULL,
                      centerx = FALSE, scalex = FALSE, prior_only = FALSE,
                      chains = 4, iter = 4e3, refresh = 500, pars = NULL,
                      control = list(adapt_delta = .9, max_treedepth = 15),
                      silent = FALSE,
                      ...) {
  if (class(family) != "family" | !family$family %in% c("binomial", "poisson")) stop ("Must provide a valid family object: binomial() or poisson().")
  if (missing(formula) | class(formula) != "formula") stop ("Must provide a valid formula object, as in y ~ x + z or y ~ 1 for intercept only.")
  if (missing(data) | missing(C)) stop("Must provide data (a data.frame or object coercible to a data.frame) and connectivity matrix C.")
  if (scalex) centerx = TRUE
  if (silent) refresh = 0
  type <- match.arg(type)
  ## GLM STUFF -------------
  a.zero <- as.array(0, dim = 1)
  tmpdf <- as.data.frame(data)
  mod.mat <- model.matrix(formula, tmpdf)
  if (nrow(mod.mat) < nrow(tmpdf)) stop("There are missing (NA) values in your data.")
  n <- nrow(mod.mat)
  ## ICAR STUFF -------------
  if (any(dim(C) != n)) stop("Dimensions of matrix C must match the number of observations. See ?shape2mat for help creating C.")
  ## GLM STUFF -------------
  family_int <- family_2_int(family)
  intercept_only <- ifelse(all(dimnames(mod.mat)[[2]] == "(Intercept)"), 1, 0)
  if (intercept_only) { # x includes slx, if any, for prior specifictions; x.list$x is the processed model matrix without slx terms.
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
            W <- C / rowSums(C)
            Wx <- SLX(f = slx, DF = tmpdf, SWM = W)
            if (scalex) Wx <- scale(Wx)            
            dwx <- ncol(Wx)
            dw_nonzero <- sum(W!=0)
            wx_idx <- as.array( which(paste0("w.", dimnames(x)[[2]]) %in% dimnames(Wx)[[2]]), dim = dwx )
            x <- cbind(Wx, x)
    }
    dbeta_prior <- ncol(x) ## dimensions of beta prior; 
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
  ## PARAMETER MODEL STUFF -------------  
  is_student <- FALSE ## family$family == "student_t"
  user_priors <- list(intercept = prior_intercept, beta = prior, sigma = NULL, nu = NULL, alpha_tau = prior_tau)
  priors <- make_priors(user_priors = user_priors, y = y, x = x, xcentered = centerx,
                        link = family$link, offset = offset)
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
    alpha_prior = priors$intercept,
    beta_prior = t(priors$beta), 
    sigma_prior = priors$sigma,
    alpha_tau_prior = priors$alpha_tau,
    t_nu_prior = priors$nu,
    family = family_int,
  ## slx data -------------    
    W = W,
    dwx = dwx,
    wx_idx = wx_idx,
  ## if TRUE, ignore data and likelihood, return prior model
    prior_only = prior_only
  )
  ## ICAR STUFF -------------
  if (inherits(scale_factor, "NULL")) inv_sqrt_scale_factor = NULL else inv_sqrt_scale_factor = 1 / sqrt(scale_factor)
  iar.list <- prep_icar_data(C, inv_sqrt_scale_factor = inv_sqrt_scale_factor)
  standata <- c(standata, iar.list)
  standata$type <- match(type, c("icar", "bym", "bym2"))  
  ## DATA MODEL STUFF -------------  
  me.list <- prep_me_data(ME, x.list$x)
  standata <- c(standata, me.list)
  ## STAN STUFF -------------    
  if (family$family == "binomial") {
      standata$y <- standata$y_int <- y[,1]
      standata$trials <- y[,1] + y[,2]
  }
  pars <- c(pars, 'intercept', 'residual', 'log_lik', 'yrep', 'fitted', 'phi', 'spatial_scale')
  if (type == "bym2") pars <- c(pars, "theta", "rho")
  if (type == "bym") pars <- c(pars, "theta", "theta_scale")
  if (standata$m) pars <- c(pars, "alpha_phi")
  if (!intercept_only) pars <- c(pars, 'beta')
  if (dwx) pars <- c(pars, 'gamma')
  if (has_re) pars <- c(pars, "alpha_re", "alpha_tau")
  if (me.list$dx_me_unbounded) pars <- c(pars, "x_true_unbounded")
  if (me.list$dx_me_bounded) pars <- c(pars, "x_true_bounded")
  priors <- priors[which(names(priors) %in% pars)]
  ## PRINT STUFF -------------    
  if (!silent) print_priors(user_priors, priors)
  ## CALL STAN -------------  
  samples <- rstan::sampling(stanmodels$icar, data = standata, iter = iter, chains = chains, refresh = refresh, pars = pars, control = control, ...)
  out <- clean_results(samples, pars, is_student, has_re, C, Wx, x.list$x, me.list$x_me_unbounded_idx, me.list$x_me_bounded_idx)
  out$data <- ModData
  out$family <- family
  out$formula <- formula
  out$slx <- slx
  out$edges <- edges(C)  
  out$re <- re_list
  out$priors <- priors
  out$scale_params <- scale_params
  if (!missing(ME)) out$ME <- ME
  out$spatial <- data.frame(par = "phi", method = toupper(type))
  class(out) <- append("geostan_fit", class(out))
  return (out)
}

