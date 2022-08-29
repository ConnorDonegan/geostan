
check_me_data <- function(ME, x.df) {    
    stopifnot(inherits(ME, "list"))
    stopifnot(inherits(ME$se, "data.frame"))
    stopifnot(inherits(ME$bounds, "numeric"))
    stopifnot(length(ME$bounds) == 2)
    stopifnot(length(ME$logit) == ncol(ME$se))
    stopifnot(inherits(ME$logit, "logical") | inherits(ME$logit, "numeric"))
    check_car_parts(ME$car_parts)
    check_me_prior(ME$prior)
    if (ME$spatial_me) stopifnot(nrow(ME$car_parts$C) == nrow(ME$se))
}

#' @noRd
check_me_prior <- function(prior) {
    stopifnot(inherits(prior, "list"))
    nms <- c("df", "location", "scale", "car_rho", "rho")
    if (!any(names(prior) %in% nms)) {
        warning("In ME$prior: dropping unused priors: ", nms[which(!names(prior) %in% nms)])
        prior <- prior[which(names(prior) %in% nms)]
    }
    if (any(unlist(lapply(prior, function(x) class(x)[1])) != "prior")) stop("ME$prior must be a list of priors (objects of class 'prior'); see ?prior for help.")
    if ("rho" %in% names(prior)) names(prior)[which(names(prior) == "rho")] <- "car_rho"
    if (!inherits(prior$df, "NULL")) stopifnot(prior$df$dist == "gamma")
    if (!inherits(prior$location, "NULL")) stopifnot(prior$location$dist == "normal")
    if (!inherits(prior$scale, "NULL")) stopifnot(prior$scale$dist == "student_t")
    if (!inherits(prior$car_rho, "NULL")) stopifnot(prior$car_rho$dist == "uniform")
}

check_car_parts <- function(car_parts) {
    if(!inherits(car_parts, "list")) stop("car_parts must be a list of data for the CAR model. See ?prep_car_data.")
    if(!all(c("Ax_w", "Ax_v", "Ax_u", "nAx_w", "Cidx", "nC", "Delta_inv", "log_det_Delta_inv", "WCAR", "lambda", "C") %in% names(car_parts))) stop("car_parts is missing at least one required part. See ?prep_car_data. Did you use cmat = TRUE and lambda = TRUE?")
    stopifnot(inherits(car_parts$C, "Matrix") | inherits(car_parts$C, "matrix"))    
}

check_sar_parts <- function(sar_parts) {
    if(!inherits(sar_parts, "list")) stop("sar_parts must be a list of data for the SAR model. See ?prep_sar_data.")
    if(!all(c("ImW_w", "ImW_v", "ImW_u", "nImW_w", "Widx", "nW", "eigenvalues_w", "n", "W", "rho_min", "rho_max") %in% names(sar_parts))) stop("sar_parts is missing at least one required part. See ?prep_sar_data.")
    stopifnot(inherits(sar_parts$W, "Matrix") | inherits(sar_parts$W, "matrix"))    
}


#' check data (x, w) dimensions and class for spatial autocorrelation measures
#' @noRd
check_sa_data <- function(x, w) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix") | inherits(w, "Matrix"))
    stopifnot(all(dim(w) == length(x)))
}

