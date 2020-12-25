#' prep_me_data
#' 
#' internal function to prepare observational error model for Stan
##' @param ME user-provided list of observational error stuff
##' @param x take design matrix from x.list$x
##' 
##' @return Partial data list to pass to Stan
##' @noRd
prep_me_data <- function(ME, x) { # for x pass in x.list$x
   # some defaults
  x.df <- as.data.frame(x)
  n <- nrow(x.df)  
  a.zero <- array(0, dim = 1)
  vec.zero <- array(0, dim = c(0, n))
  dx_me_unbounded <- 0
  dx_me_bounded <- 0
  x_me_bounded_idx = a.zero
  x_me_unbounded_idx = a.zero
  bounds <- c(0, 100)
  if (is.null(ME)) { # return items in data list ready for Stan: no ME model at all
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
          x_me_bounded = vec.zero,
          x_me_unbounded = vec.zero,
          sigma_me_bounded = vec.zero,
          sigma_me_unbounded = vec.zero,
          offset_me = rep(0, times = n),
          model_offset = 0,
          spatial_me = FALSE
      )
   return(me.list)
  }
  if (!inherits(ME, "list")) stop("ME must be a list .")
  if (length(ME$offset)) { # ME model for offset 
        offset_me <- ME$offset
        model_offset <- 1
      } else {
          offset_me <- rep(0, times = n)
          model_offset <- 0
      }
      # start building the return data list
    if (!is.null(ME$spatial)) {
        if (length(ME$spatial) != 1 | !ME$spatial %in% c(0, 1, TRUE, FALSE)) stop("ME$spatial must be logical (0, 1, TRUE, or FALSE) and of length 1.")
        spatial_me = ME$spatial
    } else {
        spatial_me <- FALSE
    }
    me.list <- list(
        offset_me = offset_me,
        model_offset = model_offset,
        spatial_me = spatial_me
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
          x_me_bounded = vec.zero,
          x_me_unbounded = vec.zero,
          sigma_me_bounded = vec.zero,
          sigma_me_unbounded = vec.zero
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
          x_me_unbounded <- as.matrix(x_me_unbounded[, x_me_unbounded_order], nrow = n)
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
 
          if (any(x_me_bounded < bounds[1]) | any(x_me_bounded > bounds[2])) stop("In ME: bounded variable has elements outside of user-provided bounds (", bounds[1], ", ", bounds[2], ")")
          # handle unused parts
          if (!dx_obs) {
              x_obs <- model.matrix(~ 0, x.df) 
              x_obs_idx <- a.zero
          }
          if (!dx_me_bounded) {
              sigma_me_bounded <- x_me_bounded <- vec.zero 
              x_me_bounded_idx <- a.zero
          }
          if (!dx_me_unbounded) {
              sigma_me_unbounded <- x_me_unbounded <- vec.zero
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
          x_me_bounded = array(t(x_me_bounded), dim = c(dx_me_bounded, n)),
          x_me_unbounded = array(t(x_me_unbounded), dim = c(dx_me_unbounded, n)),
          sigma_me_bounded = array(t(sigma_me_bounded), dim = c(dx_me_bounded, n)),
          sigma_me_unbounded = array(t(sigma_me_unbounded), dim = c(dx_me_unbounded, n))
      )
      me.list <- c(me.list, me.x.list)
 }
    return(me.list)
}


#' prep_sp_me_data
#' 
#' internal function to prepare spatial (ESF) parameters for observational error model for Stan
##' @param ME User-provided ME list
##' @param me.list output from call to prep_me_data
##' @param C Connectivity matrix. If not provided by user, it must be NA.
##' @param EV Eigenvectors. If not provided by user or previously calculated, must be NA.
##' @param x take design matrix from x.list$x
##' @param stan_esf Is this in a call to stan_esf? If so, do not return EV and dev.
##' @param silent TRUE to suppress printing of priors.
##' 
##' @return Partial data list to pass to Stan
##' @noRd
prep_sp_me_data <- function(ME, me.list, C, EV, x, stan_esf = FALSE, silent) {
  n <- nrow(x)
  if (!me.list$spatial_me) { # if no spatial data model
      slab_df_me_unbounded  <- slab_scale_me_unbounded <- scale_global_me_unbounded <- array(0, dim = 0)
      slab_df_me_bounded    <- slab_scale_me_bounded   <- scale_global_me_bounded   <- array(0, dim = 0)      
      dev = 1      
      EV <- matrix(0, nrow = n, ncol = dev)
  }
  if (me.list$spatial_me) { # if there is a spatial data model
      if ("prior_rhs" %in% names(ME)) { # if user provided prior
          rhs <- ME$prior_rhs
          if (!all(c("slab_df", "slab_scale", "scale_global", "varname") %in% names(rhs))) stop("ME$prior_rhs must be a named list of vectors 'slab_df', 'slab_scale', 'scale_global', and `varname'. List parameters in the same order as they appear in `varname.`")
          # get me_unbounded index, extract priors
          x.names <- colnames(x)
          me_ub_names <- x.names[me.list$x_me_unbounded_idx]
          me_ub_prior_idx <- match(me_ub_names, rhs$varname)
          scale_global_me_unbounded <- rhs$scale_global[me_ub_prior_idx]
          slab_df_me_unbounded <- rhs$slab_df[me_ub_prior_idx]                    
          slab_scale_me_unbounded <- rhs$slab_scale[me_ub_prior_idx]
          # get me_bounded index, extract priors
          me_b_names <- x.names[me.list$x_me_bounded_idx]
          me_b_prior_idx <- match(me_b_names, rhs$varname)
          scale_global_me_bounded <- rhs$scale_global[me_b_prior_idx]
          slab_df_me_bounded <- rhs$slab_df[me_b_prior_idx]
          slab_scale_me_bounded <- rhs$slab_scale[me_b_prior_idx]
    
          if (inherits(EV, "logical")) { # make EV if needed.
              if (!is.matrix(C)) stop("If ME$spatial is TRUE (spatial measurement error model) and you don't provide both ME$rhs_prior and EV, you need to provide connectivity matrix C.") 
              if (any(dim(C) != n)) stop("Connectivity matrix C must be n by n. See ?shape2mat for help creating C.")              
              EV <- make_EV(C)
              dev <- ncol(EV)
          } else {
              if (nrow(EV) != n) stop("data and EV have different number of rows.")              
              dev = ncol(EV)
          }
          
      } else { # else set priors here
          if (!is.matrix(C)) stop("If ME$spatial is TRUE (spatial measurement error model) and you don't provide both ME$rhs_prior and EV, you need to provide connectivity matrix C.") 
          if (any(dim(C) != n)) stop("Connectivity matrix C must be n by n. See ?shape2mat for help creating C.")
          EV.list <- make_EV(C, values = TRUE)
          EV <- EV.list$eigenvectors
          npos <- sum(EV.list$eigenvalues > 0)
          dev <- ncol(EV)
          nlinks <- length(which(C != 0))      
          E_sa <- -1 / (n-1)
          S_sa <- sqrt(2/nlinks)
          scale_ev <- sd(EV[,1])
          # calculate RHS priors for each covariate to model
          if (me.list$dx_me_unbounded) { # if there are any unbounded ME variables
              slab_df_me_unbounded <- slab_scale_me_unbounded <- scale_global_me_unbounded <- vector(mode = "numeric", length = me.list$dx_me_unbounded)
              for (i in 1:me.list$dx_me_unbounded) {
                  # slab df
                  slab_df_me_unbounded[i] <- 15
                  # slab scale
                  x.tmp <- me.list$x_me_unbounded[i,]
                  scale_x <- sd(x.tmp)
                  slab_scale_me_unbounded[i] = 0.5 * (scale_x / scale_ev)
                  # global shrinkage hyperprior
                  sa <- mc(x.tmp, C)
                  Z_sa <- (sa - E_sa) / S_sa
                  a <- (6.1808 * (Z_sa + 0.6)^0.1742)/dev^0.1298
                  b <- 3.3534/(Z_sa + 0.6)^0.1742
                  denom <- 1 + exp(2.148 - a + b)
                  p0 <- round(npos/denom)
                  if(p0 >= dev) p0 <- dev - 0.1
                  scale_global_me_unbounded[i] <- min(1, p0/(dev-p0) / sqrt(n))
                  if (!silent) {
                      message(
                          "*Setting horseshoe prior for spatial data model, variable: ",
                          colnames(x)[me.list$x_me_unbounded_idx[i]],
                          paste0(" (Moran coefficient: ", sa, ")"),                  
                          "\n global shrinkage scale:", scale_global_me_unbounded[i],
                          "\n slab scale:", slab_scale_me_unbounded[i],                  
                          "\n slab degrees of freedom:", 15, "\n"
                      )
                  }              
              }
          }
          else { # else there are no unbounded ME variables
              slab_df_me_unbounded <- slab_scale_me_unbounded <- scale_global_me_unbounded <- array(0, dim = 0)
              }             
      if (me.list$dx_me_bounded) { # if there are any bounded ME variables
          slab_df_me_bounded <- slab_scale_me_bounded <- scale_global_me_bounded <- vector(mode = "numeric", length = me.list$dx_me_bounded)
          for (i in 1:me.list$dx_me_bounded) {
              # slab df
              slab_df_me_bounded[i] <- 15
              x.tmp <- me.list$x_me_bounded[i,]
              # slab scale
              scale_x <- sd(x.tmp)
              slab_scale_me_bounded[i] = 0.5 * (scale_x / scale_ev)
              # global shrinkage hyperprior
              sa <- mc(x.tmp, C)
              Z_sa <- (sa - E_sa)/S_sa
              a <- (6.1808 * (Z_sa + 0.6)^0.1742)/dev^0.1298
              b <- 3.3534/(Z_sa + 0.6)^0.1742
              denom <- 1 + exp(2.148 - a + b)
              p0 <- round(npos/denom)
              if(p0 >= dev) p0 <- dev - 0.1
              scale_global_me_bounded[i] <- min(1, p0/(dev-p0) / sqrt(n))
              if (!silent) {
                  message(
                      "*Setting horseshoe prior for spatial data model, variable: ",
                      colnames(x)[me.list$x_me_bounded_idx[i]],
                      paste0(" (Moran coefficient: ", sa, ")"),
                      "\n Global shrinkage (scale_global): ", round(scale_global_me_bounded[i], 4),
                      "\n Slab scale: ", round(slab_scale_me_bounded[i], 4),                 
                      "\n Slab degrees of freedom: ", 15, "\n"
                  )
              }          
          }
      } else { # else there are no bounded ME variables
          slab_df_me_bounded <- slab_scale_me_bounded <- scale_global_me_bounded <- array(0, dim = 0)
      }      
      }
  }
  sp.me <- list(
      slab_df_me_unbounded = array(slab_df_me_unbounded),
      slab_df_me_bounded = array(slab_df_me_bounded),
      slab_scale_me_unbounded = array(slab_scale_me_unbounded),
      slab_scale_me_bounded = array(slab_scale_me_bounded),
      scale_global_me_unbounded = array(scale_global_me_unbounded),
      scale_global_me_bounded = array(scale_global_me_bounded)
  )
  if (!stan_esf) {
      sp.me$dev <- dev
      sp.me$EV <- EV
  }
  return(sp.me)
}

