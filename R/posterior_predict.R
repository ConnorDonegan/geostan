#' Sample from the posterior predictive distribution
#'
#' @description Draw samples from the posterior predictive distribution of a fitted \code{geostan} model. Use the original data or new data, with or without the spatial component and random effects terms.
#' @export
#' @param object A \code{geostan_fit} object.
#' @param newdata A \code{data.frame} or \code{matrix} with data to use for predictions. Must be named and the names must match the model formula. 
#' @param W Spatial connectivity matrix for calculating spatially lagged covariates (slx), if needed. If the model has an slx term then \code{newdata} must have the same number of rows as the data used to fit the model.
#' @param samples Number of samples to take from the posterior distribution. Maximum is the total number of samples stored in the model.
#' @param predictive Return the predictive distribution? Defaults to \code{TRUE}. If \code{FALSE}, then the linear predictor is returned after applying the inverse link function. 
#' @param re_form If \code{re_form = NA} any random effects terms will be ignored when making predictions. If default \code{re_form = NULL} is set, then the results will include random effects terms. Those are the only options; this syntax is used to maintain some consistency with other R packages which have additional options.
#' @param spatial Include the spatial component in the model predictions? Defaults to \code{TRUE} and will be ignored if the model does not have a spatial comonent (i.e. \code{\link[geostan]{stan_glm}}). For models fit by \code{\link[geostan]{stan_esf}}, \code{\link[geostan]{stan_icar}}, and \code{\link[geostan]{stan_bym2}} this option requires that \code{newdata} have the same number of observations as the data that the model was fit to.
#' @param seed A single integer value to be used in a call to \code{\link[base]{set.seed}} before taking samples from the posterior distribution. Passing a value to \code{seed} enables you to obtain the same results each time (by using the same seed).
#' @param centerx Should \code{newdata} be centered using the means of the variables in the original data used to fit the model (stored in \code{fit$scale_params})? Defaults to \code{TRUE} but is (effectively) ignored if the model was fit with \code{centerx = FALSE}.
#' @param scalex Should \code{newdata} be scaled using the standard deviations of the variables in the original data used to fit the model (stored in \code{fit$scale_params})? Defaults to \code{FALSE}.
#' 
#' @return A matrix of size \code{S} x \code{N} containing samples from the posterior predictive distribution, where \code{S} is the number of samples and \code{N} is the number of observations. It is of class \code{matrix} and \code{ppd}. 
#'
#' @examples
#' library(bayesplot)
#' library(ggplot2)
#' library(sf)
#' data(ohio)
#' ## fit Ohio election model with a spatial filter
#' fit <- stan_esf(gop_growth ~ log(pop_density) + historic_gop + college_educated,
#'                 data = ohio,
#'                 family = student_t(),
#'                 C = shape2mat(ohio, "B"),
#'                 iter = 400, chains = 1)
#' 
#' ## posterior predictive check (density overlay)
#' ## compare distribution of observations to the predictive distribution
#' y <- ohio$gop_growth
#' yrep <- posterior_predict(fit, samples = 150)
#' ppc_dens_overlay(y, yrep)
#'
#' ## model convict-leasing era sentencing risk in Florida
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' fit <- stan_esf(sents ~ offset(expected_sents),
#'                 re = ~ name,
#'                 data = sentencing,
#'                 family = poisson(),
#'                 C = C,
#'                 chains = 1, iter = 400)
#' ## posterior predictive checks
#' yrep <- posterior_predict(fit, samples = 200, seed = 1)
#'  # density overlay
#' ppc_dens_overlay(sentencing$sents, yrep)
#'  # map realizations from the model, as relative risk  
#' fl <- st_as_sf(sentencing)
#' E <- fl$expected_sents
#' map_pp <- function(i) {
#'    ggplot(fl) +
#'    geom_sf(aes(fill = yrep[i,] / E) +
#'    scale_fill_gradient2(midpoint = 1) +
#'    theme_void()
#' }
#' map_pp(1)
#' map_pp(20)
#' map_pp(90)
#'
posterior_predict <- function(object, newdata, W, samples, predictive = TRUE, re_form = NULL, spatial = TRUE, seed, centerx = TRUE, scalex = FALSE) {
    if (!inherits(object, "geostan_fit")) stop ("object must be of class geostan_fit.")
    N <- nrow(as.matrix(object))                                 
    if (!missing(seed)) set.seed(seed)                                      
    if (missing(samples)) samples <- N
    if (samples > N) {
        warning (paste0("Cannot draw more samples than were taken from the posterior. Using samples = ", N))
        samples <- N
    }
    idx <- sample(N, samples)
    if (missing(newdata)) {
        preds <- as.matrix(object, pars = "yrep")[idx,]
        class(preds) <- append("ppd", class(preds))
        return(preds)
    } else {
        newdata <- as.data.frame(newdata)
    }
    if (class(object$slx) == "formula") {
        if (missing(W)) stop ("You must provide spatial weights matrix W to calculate spatial lag of X terms (slx)")
        if ( nrow(newdata) != nrow(W) | (nrow(newdata) != nrow(object$data)) ) stop ("W, newdata, and object$data must have same number of rows.")
    }
    family <- object$family$family
    link <- object$family$link
    if (family == "binomial") {
        y <- model.response(model.frame(object$formula, newdata))
        trials <- as.integer(y[, 1] + y[, 2])
    } else {
        y <- as.character(object$formula[[2]])
        if (!y %in% names(newdata))  newdata[,y] <- 0  
    }
    x <- model.matrix(object$formula, data = newdata) 
    x <- remove_intercept(x)
    if (centerx) centerx <- fit$scale_params$center
    if (scalex) scalex <- fit$scale_params$scale
    x <- scale(x, center = centerx, scale = scalex)
    offset <- model.offset(model.frame(object$formula, newdata))
    if (!is.null(offset) & family == "poisson") offset <- log(offset)
    if ( spatial & (class(object$slx) == "formula") ) {
         Wx <- SLX(f = object$slx, DF = newdata, SWM = W, cx = FALSE, sx = FALSE)
         x <- cbind(Wx, x)
    }
    alpha <- as.matrix(object, pars = "intercept")[idx,]
    beta <- as.matrix(object, pars = "beta")[idx,]
    mu <- alpha + beta %*% t(x)
    if (!is.null(offset)) mu <- sweep(mu, 2, offset, "+")
    if (is.null(re_form) & !is.na(object$re[1])) {
        newdata$order <- 1:nrow(newdata)
        re_term <- as.character( object$re$formula[[2]] )
        re_idx <- which( names(newdata) == re_term )
        if ( !all(newdata[,re_idx] %in% object$re$data$id) ) stop ("New levels for the random effects term are not allowed. Find fitted levels in object$re$data$id.")
        names(newdata)[re_idx] <- "id"
        newdata <- merge(newdata, object$re$data, by = "id", sort = FALSE)
        newdata <- newdata[order(newdata$order),]
        alpha_re <- as.matrix(object, pars = "alpha_re")[idx,]
        alpha_re <- alpha_re[,paste0("alpha_re[",newdata$idx,"]")]
        mu <- mu + alpha_re
    }
    if (spatial & !object$spatial$method %in% c("None", "none", "Exchangeable")) { 
       if (nrow(object$data) != nrow(newdata)) stop ("If spatial = TRUE, newdata must contain the same number of observations as the original data (i.e. nrow(object$data)).")
       sp <- spatial(object, summary = FALSE)[idx,]
       mu <- mu + sp
    }
    if (link == "log") mu <- exp(mu) 
    if (link == "logit") mu <- inv_logit(mu)
    if (!predictive) return (mu)
    if (family == "gaussian") {
        sigma <- as.matrix(object, pars = "sigma")[idx,]
        preds <- .pp_gaussian(mu, sigma)
    }
    if (family == "student_t") {
        sigma <- as.matrix(object, pars = "sigma")[idx,]
        nu <- as.matrix(object, pars = "nu")[idx,]
        preds <- .pp_student(nu, mu, sigma)
    }
    if (family == "poisson") preds <- .pp_poisson(mu)
    if (family == "binomial") preds <- .pp_binomial(mu, trials)
    class(preds) <- append("ppd", class(preds))
    return(preds)
}


inv_logit <- function(x) exp(x)/(1+exp(x))

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}

.pp_student <- function(nu, mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
   mu[s,] + rt(ncol(mu), df = nu[s]) * sigma[s]
  }))
}

.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s, ])
  }))
}

.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s, ]) 
  }))
}

