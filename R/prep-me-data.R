#' internal function to prepare observational error data for Stan
##' @param ME user-provided list of observational error stuff
##' @param family user-provided family object
##' @param x take design matrix from x.list$x
##' 
##' @return Partial data list to pass to Stan
##' @noRd
prep_me_data <- function(ME, family, x) { # for x pass in x.list$x
   # some defaults
  a.zero <- as.array(0, dim = 1)
  dx_me_unbounded <- 0
  dx_me_bounded <- 0
  x_me_bounded_idx = a.zero
  x_me_unbounded_idx = a.zero
  bounds <- c(0, 100)
  x.df <- as.data.frame(x)
  n <- nrow(x.df)  
      # return items in data list ready for Stan: no ME model at all
  if (is.null(ME)) {
      x_obs <- x
      dx_obs <- ncol(x_obs)
      if (dx_obs) {
          x_obs_idx <- as.array(1:dx_obs, dim = dx_obs)
      } else {
          x_obs_idx <- a.zero
      }
      me.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,          
          x_obs = x_obs,
          x_me_bounded = matrix(0, nrow = n, ncol = 1),
          x_me_unbounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_bounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_unbounded = matrix(0, nrow = n, ncol = 1),
          offset_me = rep(0, times = n),
          model_offset = 0
      )
   return(me.list)
  }
  if (!inherits(ME, "list")) stop("ME must be a list .")
  if (length(ME$offset)) { # ME model for offset #
        offset_me <- ME$offset
        model_offset <- 1
      } else {
          offset_me <- rep(0, times = n)
          model_offset <- 0
      }
      # start building the return data list
    me.list <- list(
        offset_me = offset_me,
        model_offset = model_offset
    )
  if (is.null(ME$se)) {
      # in this case there is just an offset me model but not ME.X
      x_obs <- x
      dx_obs <- ncol(x_obs)
          if (dx_obs) {
              x_obs_idx <- as.array(1:dx_obs, dim = dx_obs)
          } else {
              x_obs_idx <- a.zero
          }
      me.x.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,          
          x_obs = x_obs,
          x_me_bounded = matrix(0, nrow = n, ncol = 1),
          x_me_unbounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_bounded = matrix(0, nrow = n, ncol = 1),
          sigma_me_unbounded = matrix(0, nrow = n, ncol = 1)
      )
          me.list <- c(me.list, me.x.list)
 } else {
     ## now deal with observational error in covariates, if any
        if (!inherits(ME$se, "data.frame")) stop("ME must be a list in which the element named ME is of class data.frame, containing standard errors for the observations.")
        if  (!all(names(ME$se) %in% names(x.df))) stop("All column names in ME$se must be found in the model matrix (from model.matrix(formula, data)). This error may occur if you've included some kind of data transformation in your model formula, such as a logarithm or polynomial, which is not supported for variables with sampling/measurement error.")
        if (length(ME$bounded)) {
            if (length(ME$bounded) != ncol(ME$se)) stop("ME mis-specified: bounded must be a vector with one element per column in the ME dataframe.")
            bounded <- which(ME$bounded == 1)
            not.bounded <- which(ME$bounded == 0)
            if (length(ME$bounds)) {
                if(length(ME$bounds) != 2 | !inherits(ME$bounds, "numeric")) stop("ME$bounds must be numeric vector of length 2.")
                bounds <- ME$bounds
            }
        } else {
            bounded <- integer(0) #rep(0, times = ncol(ME$se))
            not.bounded <- 1:ncol(ME$se) #rep(1, times = ncol(ME$se))
        }           
           # gather any/all variables without ME
          x_obs_idx <- as.array( which( !names(x.df) %in% names(ME$se) )) 
          x_obs <- as.data.frame(x.df[, x_obs_idx])
          dx_obs <- ncol(x_obs)
          # now X.me needs to be parsed into bounded/non-bounded variables and ordered as x
          nm_me_unbounded <- names(ME$se)[not.bounded]
          x_me_unbounded <- data.frame( x.df[, nm_me_unbounded] )
          names(x_me_unbounded) <- nm_me_unbounded
          x_me_unbounded_order <- na.omit( match(names(x.df), names(x_me_unbounded)) )
          x_me_unbounded <- data.frame(x_me_unbounded[, x_me_unbounded_order])
          dx_me_unbounded <- ncol(x_me_unbounded)
          sigma_me_unbounded <- as.matrix(ME$se[,not.bounded], nrow = n)
          sigma_me_unbounded <- as.matrix(sigma_me_unbounded[, x_me_unbounded_order], nrow = n)
          x_me_unbounded_idx <- as.array( which( names(x.df) %in% nm_me_unbounded ))
          
          nm_me_bounded <- names(ME$se)[bounded]
          x_me_bounded <- data.frame( x.df[, nm_me_bounded] )
          names(x_me_bounded) <- nm_me_bounded
          x_me_bounded_order <- na.omit( match(names(x.df), names(x_me_bounded)) )
          x_me_bounded <- as.matrix(x_me_bounded[, x_me_bounded_order], nrow = n)
          dx_me_bounded <- ncol(x_me_bounded)
          sigma_me_bounded <- as.matrix(ME$se[,bounded], nrow = n)
          sigma_me_bounded <- as.matrix(sigma_me_bounded[, x_me_bounded_order], nrow = n)
          x_me_bounded_idx <- as.array( which( names(x.df) %in% nm_me_bounded ))
          # handle unused parts
          if (!dx_obs) {
              x_obs <- model.matrix(~ 0, x.df) 
              x_obs_idx <- a.zero
          }
          if (!dx_me_bounded) {
              sigma_me_bounded <- x_me_bounded <- matrix(0, nrow = n, ncol = 1)
              x_me_bounded_idx <- a.zero
          }
          if (!dx_me_unbounded) {
              sigma_me_unbounded <- x_me_unbounded <- matrix(0, nrow = n, ncol = 1)
              x_me_unbounded_idx <- a.zero
          }
                      # return items in data list ready for Stan: with ME model for covariates
          me.x.list <- list(
          dx_obs = dx_obs,
          dx_me_unbounded = dx_me_unbounded,
          dx_me_bounded = dx_me_bounded,
          x_obs_idx = x_obs_idx,
          x_me_bounded_idx = x_me_bounded_idx,
          x_me_unbounded_idx = x_me_unbounded_idx,
          bounds = bounds,
          x_obs = x_obs,
          x_me_bounded = x_me_bounded,
          x_me_unbounded = x_me_unbounded,
          sigma_me_bounded = sigma_me_bounded,
          sigma_me_unbounded = sigma_me_unbounded
      )
      me.list <- c(me.list, me.x.list)
      }   
    return(me.list)
}
