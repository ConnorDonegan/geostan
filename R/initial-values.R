#' Construct a function that will draw initial values for ME models
#' 
#' @noRd
#'
#' @param FRAME_NUMBER the sys.nframe() indexing the environment immediately within the stan_* function call
#' 
#' @description When logit transform is used with ME models, the initial values need to remain within zero-one range. This will do so safely regardless of user-provided parameter constraints (bounds), enabling greater flexibility in the types of covariates that can be included in the same model.
init_fn_builder <- function(FRAME_NUMBER) {
    logits = sys.frame(which=FRAME_NUMBER)$standata$use_logit
    lims = sys.frame(which=FRAME_NUMBER)$standata$bounds                
    nvar = sys.frame(which=FRAME_NUMBER)$standata$dx_me
    nobs = sys.frame(which=FRAME_NUMBER)$standata$n
    inits_fn <- function(
                         bounds = lims,
                         K = nvar,
                         N = nobs,
                         use_logit = logits) {
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
        l <- list(
            x_true = x_true_init
        )
        return (l)
    }
    return (inits_fn)
}
