


#' @description Return index of observed, censored y; elsewhere, use results to replace NAs with zeros
#' This will stop if there are missing values and the censor_point argument is not being used; outside of this call, must check that censor_point argument is only used with Poisson likelihood.
#' 
#' @param censor the censor_point argument
#' @param frame from model.frame(formula, tmpdf, na.action = NULL)
#' @noRd
#' 
handle_censored_y <- function(censor, frame) {
    y_raw <- model.response(frame)
    if (inherits(y_raw, "matrix")) {
        y_raw <- y_raw[,1] + y_raw[,2]
    }
    y_mis_idx <- which(is.na(y_raw))
    y_obs_idx <- which(!is.na(y_raw))
    return (list(n_mis = length(y_mis_idx),
                n_obs = length(y_obs_idx),
                y_mis_idx = y_mis_idx,
                y_obs_idx = y_obs_idx,
                censor_point = censor)
           )
}

#' @noRd
handle_missing_x <- function(frame) {
    check_NAs <- apply(as.matrix(frame[,-1]), 2, function(c) any(is.na(c)))
    if (any(check_NAs)) stop("Missing values found in covariates or offset term. Missing values are only permitted in the outcome variable.")
}
