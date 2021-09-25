

#' Prior distributions
#' 
#' @details
#' 
#' The prior distribution functions are used to set the values of prior parameters.
#'
#' Users can control the values of the parameters, but the distribution (model) itself is fixed. The intercept and regression coefficients are given Gaussian prior distributions and scale parameters are assigned Student's t prior distributions. Degrees of freedom pararmeters are assigned gamma priors, and the spatial autocorrelation parameter in the CAR model, rho, is assigned a uniform prior. The horseshoe (`hs`) model is used by \code{\link[geostan]{stan_esf}}.
#'
#' Note that the `variable` argument is used internally by `geostan`, and any user provided values will be ignored.
#'
#' ### Parameterizations
#'
#' For details on how any distribution is parameterized, see the Stan Language Functions Reference document: \url{https://mc-stan.org/users/documentation/}.
#'
#' @return An object of class `prior` which will be used internally by **geostan** to set parameters of prior distributions. 
#' 
#' @examples
#'
#' prior <- list()
#' prior$beta <- normal(c(0, 0), c(1, 1))
#' prior$intercept <- normal(-5, 3)
#' \dontrun{
#' fit <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + ICE + college,
#'                 re = ~ GEOID,
#'                 data = georgia,
#'                 family = poisson(),
#'                 prior = prior,
#'                 prior_only = TRUE)
#' plot(fit)
#' }
#'
#' ME <- list()
#' ME$se <- data.frame(insurance = georgia$insurance.se)
#' ME$prior <- list()
#' ME$prior$nu <- gamma(3, 0.2)
#' ME$prior$location <- normal(50, 50)
#' ME$prior$scale <- student_t(12, 10, 20)
#' \dontrun{
#' fit <- stan_glm(log(rate.male) ~ insurance, 
#'                 data = georgia,
#'                 ME = ME,
#'                 prior_only = TRUE)
#' }
#' @name priors
#' @md
#' 
NULL

#' @rdname priors
#' @param lower,upper lower and upper bounds of the distribution
#' @param variable A reserved slot for the variable name; if provided by the user, this may be ignored by **geostan**. 
#' @export
#' @md
uniform <- function(lower, upper, variable = NULL) {
    out <- list(dist = "uniform", lower = lower, upper = upper, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)    
}


#' @param location Location parameter(s), numeric value(s)
#' @param scale  Scale parameter(s), positive numeric value(s)
#' 
#' @rdname priors
#' @export
#' @md
normal <- function(location = 0, scale, variable = NULL) {
  validate_positive_parameter(scale)
  out <- list(dist = "normal", location = location, scale = scale, variable = variable)
  class(out) <- append("prior", class(out))
  return (out)
}

#'
#' @param df Degrees of freedom, positive numeric value(s)
#' 
#' @return
#'
#' ### Student's t
#' 
#' Return value for \code{student_t} depends on the input; if no arguments are provided (specifically, if the scale parameter is missing), this will return an object of class 'family'; if at least the scale parameter is provided, `student_t` will return an object of class `prior` containing parameter values for the Student's t distribution.
#' 
#' @examples
#' \dontrun{
#' fit = stan_glm(log(rate.male) ~ 1, data = georgia, family = student_t())
#' }
#' @rdname priors
#' @export
#' @md
student_t <- function(df = 10, location = 0, scale, variable = NULL) {
    if (missing(scale)) {
        family <- list(family = "student_t", link = 'identity')
        class(family) <- "family"
        return(family)
    }
    validate_positive_parameter(scale)
    validate_positive_parameter(df)
    out <- list(dist = "student_t", df = df, location = location, scale = scale, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)
}

#' @param alpha shape parameter, positive numeric value(s)
#' @param beta inverse scale parameter, positive numeric value(s)
#'
#' 
#' @rdname priors
#' @export
#' @md
gamma <- function(alpha, beta, variable = NULL) {
    out <- list(dist = "gamma", alpha = alpha, beta = beta, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)    
}


#' @param global_scale Control the (prior) degree of sparsity in the horseshoe model (0 < global_scale < 1).
#' @param slab_df Degrees of freedom for the Student's t model for large coefficients in the horseshoe model (slab_df > 0).
#' @param slab_scale Scale parameter for the Student's t model for large coefficients in the horseshoe model (slab_scale > 0).
#' @details
#'
#' ### The horseshoe prior
#' 
#' The horseshoe prior is used by \code{\link[geostan]{stan_esf}} as a prior for the eigenvector coefficients. The horseshoe model encodes a prior state of knowledge that effectively states, 'I believe a small number of these variables may be important, but I don't know which of them is important.' The horseshoe is a normal distribution with unkown scale (Polson and Scott 2010):
#' ```
#'        beta_j ~ Normal(0, tau^2 lambda_j^2)
#' ```
#' The scale parameter for this prior is the product of two terms: `lambda_j^2` is specific to the variable `beta_j`, and `tau^2` is known as the global shrinkage parameter.
#'
#' The global shrinkage parameter is assigned a half-Cauchy prior:
#' ```
#'        tau ~ Cauchy(0, global_scale * sigma)
#' ```
#' where `global_scale` is provided by the user and `sigma` is the scale parameter for the outcome variable; for Poisson and binomial models, sigma is fixed at one. Use `global_scale` to control the overall sparsity of the model.
#'
#' The second part of the model is a Student's t prior for `lambda_j`. Most `lambda_j` will be small, since the model is half-Cauchy:
#' ```
#'        lambda_j ~ Cauchy(0, 1)
#' ```
#' This model results in most `lambda_j` being small, but due to the long tails of the Cauchy distribution, strong evidence in the data can force any particular `lambda_j` to be large. Piironen and Vehtari (2017) adjust the model so that those large `lambda_j` are effectively assigned a Student's t model:
#' ```
#'        Big_lambda_j ~ Student_t(slab_df, 0, slab_scale)
#' ```
#' This is a schematic representation of the model; see Piironen and Vehtari (2017) or Donegan et al. (2020) for details.
#' 
#' @source
#'
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. *Spatial Statistics*. \doi{10.1016/j.spasta.2020.100450} (open access: \doi{10.31219/osf.io/fah3z}).
#' 
#' Polson, N.G. and J.G. Scott (2010). Shrink globally, act locally: Sparse Bayesian regularization and prediction. *Bayesian Statistics* 9, 501-538. \doi{10.1093/acprof:oso/9780199694587.003.0017}.
#'
#' Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In *Electronic Journal of Statistics*, 11(2):5018-5051. \doi{10.1214/17-EJS1337SI}.
#' 
#' @rdname priors
#' @export
#' @md
hs <- function(global_scale = 1, slab_df = 10, slab_scale, variable = "beta_ev") {
    validate_positive_parameter(global_scale)
    validate_positive_parameter(slab_df)
    validate_positive_parameter(slab_scale)
    out <- list(dist = "hs", global_scale = global_scale, slab_df = slab_df, slab_scale = slab_scale, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)    
}

#' @noRd
#' @param x The value to check.
#' @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_positive_parameter <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}


#' @noRd
#' 
append_priors <- function(standata, priors_made) {
    dl <- list(
        prior_alpha = c(priors_made$intercept$location, priors_made$intercept$scale),
        dbeta_prior = length(priors_made$beta$location),
        prior_beta_location = as.array(priors_made$beta$location),
        prior_beta_scale = as.array(priors_made$beta$scale),
        prior_sigma = c(priors_made$sigma$df, priors_made$sigma$location, priors_made$sigma$scale),
        prior_alpha_tau = c(priors_made$alpha_tau$df, priors_made$alpha_tau$location, priors_made$alpha_tau$scale),
        prior_t_nu = c(priors_made$nu$alpha, priors_made$nu$beta)
    )
    return(c(standata, dl))
}


#' @export
#' @noRd
#' @param x a prior distribution object
#' @param digits number of digits to print
#' @param ... additional arguments passed to \code{\link[base]{print.data.frame}}
#' @method print prior
print.prior <- function(x, digits = 2, ...) {
    nm <- x$dist
    var <- ifelse(is.null(x$variable), FALSE, x$variable)
    message("Distribution: ", nm)
    if (nm == "gamma")     df <- as.data.frame(x[c('alpha', 'beta')])
    if (nm == "student_t") df <- as.data.frame(x[c('df', 'location', 'scale')])
    if (nm == "uniform")   df <- as.data.frame(x[c('lower', 'upper')])
    if (nm == "normal")    df <- as.data.frame(x[c('location', 'scale')])
    if (nm == "hs")        df <- as.data.frame(x[c('global_scale', 'slab_df', 'slab_scale')])
    print(df, digits = digits, ...)           
}
