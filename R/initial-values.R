#' Construct a function that will draw initial values for ME models
#' 
#' @noRd
#'
#' @param FRAME_NUMBER the sys.nframe() indexing the environment immediately within the stan_* function call
#' 
#' @description When logit transform is used with ME models, the initial values need to remain within zero-one range. This will do so safely regardless of user-provided parameter constraints (bounds), enabling greater flexibility in the types of covariates that can be included in the same model.
#'
#' The adjustment for censored observations is to keep the intercept (mean log-rate) below zero (centered on the mean value, which should not evaluate to log(0)), and if its  CAR model, to also initiate log-rates at less than zero. Censored models are only available for Poisson models.
#'
#' @importFrom truncnorm rtruncnorm
init_fn_builder <- function(FRAME_NUMBER) {
    logits = sys.frame(which=FRAME_NUMBER)$standata$use_logit
    lims = sys.frame(which=FRAME_NUMBER)$standata$bounds                
    nvar = sys.frame(which=FRAME_NUMBER)$standata$dx_me
    nobs = sys.frame(which=FRAME_NUMBER)$standata$n
    censor_point = sys.frame(which=FRAME_NUMBER)$standata$censor_point
    alpha_prior = sys.frame(which=FRAME_NUMBER)$standata$prior_alpha
    is_car = as.logical(sys.frame(which=FRAME_NUMBER)$standata$car)
    is_poisson = sys.frame(which=FRAME_NUMBER)$standata$family == 3
    is_car_poisson = as.logical(is_car & is_poisson)
    inits_fn <- function(
                         bounds = lims,
                         K = nvar,
                         N = nobs,
                         use_logit = logits,
                         censor = censor_point > 0,
                         prior_mean = alpha_prior[1],
                         start_log_lambda = is_car_poisson 
                         ) {
        l <- NULL
        if (any(use_logit == 1)) {
            S <- matrix(0, nrow = K, ncol = N)        
            for (j in 1:K) {
            if (use_logit[j]) {
                lwr = 0; upr = 1;
            } else {
                lwr = bounds[1]; upr = bounds[2];
            }
            S[j,] <- truncnorm::rtruncnorm(n = N,
                                           a = lwr,
                                           b = upr,
                                           mean = 0,
                                           sd = 1)        
            }
            x_true_init <- array(S, dim = c(K, N))
            l1 <- list(
                x_true = x_true_init
            )
            l <- c(l, l1)
        }
        if (censor > 0) {
            a <- rnorm(n = 1, mean = prior_mean, sd = 0.25)
            l2 <- list(intercept = a)
            l <- c(l, l2)
            if (start_log_lambda) {
                phi <- rnorm(n = N, mean = prior_mean, sd = 0.25)
                tau <- truncnorm::rtruncnorm(n = 1, a = 0, mean = 0, sd = 0.25)
                l3 <- list(log_lambda = phi, log_lambda_mu = phi, car_scale = array(tau))
                l <- c(l, l3)
             }
        }        
        return (l)
    }
    return (inits_fn)
}
