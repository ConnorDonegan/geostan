

#' Prior distributions
#' @name priors
#' @details
#' The prior distribution functions are used to set the values of prior parameters.
#'
#' Note that users can control the values of the parameters, but the distribution (model) itself is fixed. The intercept and regression coefficients are given Gaussian prior distributions and scale parameters are assigned Student's t prior distributions. Degrees of freedom pararmeters are assigned gamma priors, and the spatial autocorrelation parameter in the CAR model, rho, is assigned a uniform prior. The horseshoe (`hs`) model is used by \code{\link[geostan]{stan_esf}}.
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
NULL

#' @rdname priors
#' @export
uniform <- function(lower, upper, variable = NULL) {
    out <- list(dist = "uniform", lower = lower, upper = upper, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)    
}

#' @rdname priors
#' @export
#'
#' @details
#'
#' Note that the `variable` argument is used internally by `geostan`, and any user provided values will be ignored.
#' 
normal <- function(location = 0, scale, variable = NULL) {
  validate_positive_parameter(scale)
  out <- list(dist = "normal", location = location, scale = scale, variable = variable)
  class(out) <- append("prior", class(out))
  return (out)
}

#' @rdname priors
#' @return Return value for \code{student_t} depends on the input; if no arguments are provided (specifically, if the scale parameter is missing), this will return an object of class 'family'; if at least the scale parameter is provided, `student_t` will return an object of class `prior` containing parameter values for the Student's t distribution.
#' @export
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

#' @rdname priors
#' @export
gamma <- function(alpha, beta, variable = NULL) {
    out <- list(dist = "gamma", alpha = alpha, beta = beta, variable = variable)
    class(out) <- append("prior", class(out))
    return (out)    
}

#' @rdname priors
#' @export
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
#' @method print prior
print.prior <- function(prior) {
    nm <- prior$dist
    var <- ifelse(is.null(prior$variable), FALSE, prior$variable)
    message("Distribution: ", nm)
    if (nm == "gamma")     df <- as.data.frame(prior[c('alpha', 'beta')])
    if (nm == "student_t") df <- as.data.frame(prior[c('df', 'location', 'scale')])
    if (nm == "uniform")   df <- as.data.frame(prior[c('lower', 'upper')])
    if (nm == "normal")    df <- as.data.frame(prior[c('location', 'scale')])
    if (nm == "hs")        df <- as.data.frame(prior[c('global_scale', 'slab_df', 'slab_scale')])
    print(df)           
}
