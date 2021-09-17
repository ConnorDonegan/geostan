
check_me_data <- function(ME, x.df) {    
    stopifnot(inherits(ME, "list"))
    stopifnot(inherits(ME$se, "data.frame"))
    if (!all(names(ME$se) %in% names(x.df))) stop("All column names in ME$se must be found in the model matrix (from model.matrix(formula, data)). This error may occur if you've included some kind of data transformation in your model formula, such as a logarithm or polynomial.")
}

check_car_parts <- function(car_parts) {
    if(!inherits(car_parts, "list")) stop("car_parts must be a list of data for the CAR model. See ?prep_car_data.")
    if(!all(c("Ax_w", "Ax_v", "Ax_u", "nAx_w", "Cidx", "nC", "Delta_inv", "log_det_Delta_inv", "WCAR", "lambda", "C") %in% names(car_parts))) stop("car_parts is missing at least one required part. See ?prep_car_data. Did you use cmat = TRUE and lambda = TRUE?")
    stopifnot(inherits(car_parts$C, "Matrix") | inherits(car_parts$C, "matrix"))    
}

#' check data (x, w) dimensions and class for spatial autocorrelation measures
#' @noRd
check_sa_data <- function(x, w) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix") | inherits(w, "Matrix"))
    stopifnot(all(dim(w) == length(x)))
}
