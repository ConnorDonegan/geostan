#' prep_me_data
#' 
#' internal function to prepare observational error model for Stan
#' 
##' @param ME user-provided list of observational error stuff
##' 
##' @param x design matrix, excluding any Wx (slx) terms
##' 
##' @return Partial data list to pass to Stan. 
##' @noRd
prep_me_data <- function(ME, x) { # for x pass in x_no_Wx !
    x.df <- as.data.frame(x)
    n <- nrow(x.df)  
    dx_me <- 0
    x_me_idx = a.zero()
    bounds <- c(-Inf, Inf)
    if (is.null(ME)) { # return items in data list ready for Stan: no ME model at all. 
        x_obs <- x
        dx_obs <- ncol(x_obs)
        if (dx_obs) {
            x_obs_idx <- as.array(1:dx_obs, dim = dx_obs)
        } else {
            x_obs_idx <- a.zero()
        }
        me.list <- list(
            dx_obs = dx_obs,
            dx_me = dx_me,
            x_obs_idx = x_obs_idx,
            x_me_idx = x_me_idx,
            bounds = bounds,          
            x_obs = x_obs,
            x_me = a.n.zeros(n),
            sigma_me = a.n.zeros(n),
            ME_prior_car_rho = vec.n.zeros(2),
            has_me = FALSE,
            spatial_me = FALSE          
        )
        me.list <- c(me.list, car_parts_shell(n))
        pl <- me_priors(ME, me.list)
        me.list <- c(me.list, pl)
        return(me.list)
    } #/#
    check_me_data(ME, x.df)    
    if (length(ME$bounds)) {
        if(length(ME$bounds) != 2 | !inherits(ME$bounds, "numeric")) stop("ME$bounds must be numeric vector of length 2.")
        bounds <- ME$bounds        
    } else {
        bounds <- c(-Inf, Inf)
    }           
    # gather any/all variables without ME
    x_obs_idx <- as.array( which( !names(x.df) %in% names(ME$se) )) 
    x_obs <- as.data.frame(x.df[, x_obs_idx])
    dx_obs <- ncol(x_obs)
    # gather ME variables
    x_me_idx <- as.array(na.omit(match(names(ME$se), names(x.df))))
    x_me <- as.data.frame(x.df[,x_me_idx])
    dx_me <- ncol(x_me)
    sigma_me <- ME$se    
    if (any(x_me < bounds[1]) || any(x_me > bounds[2])) stop("In ME: bounded variable has elements outside of user-provided bounds (", bounds[1], ", ", bounds[2], ")")
    # handle unused parts
    if (!dx_obs) {
        x_obs <- model.matrix(~ 0, x.df) 
        x_obs_idx <- a.zero()
    }
    # return items in data list ready for Stan: with ME model for covariates
    me.list <- list(
        dx_obs = dx_obs,
        dx_me = dx_me,
        x_obs_idx = x_obs_idx,
        x_me_idx = x_me_idx,
        bounds = bounds,
        x_obs = x_obs, 
        x_me = array(t(x_me), dim = c(dx_me, n)),
        sigma_me = array(t(sigma_me), dim = c(dx_me, n)),
        nm_me = names(ME$se),
        has_me = TRUE
    )
    if (inherits(ME$car_parts, "NULL")) {
        me.list$spatial_me <- FALSE        
        me.list <- c(me.list, car_parts_shell(n))
    } else {        
        check_car_parts(ME$car_parts)
        me.list$spatial_me <- TRUE
        me.list <- c(me.list, ME$car_parts)
    }
    pl <- me_priors(ME, me.list)
    me.list <- c(me.list, pl)
    return(me.list)
}


me_priors <- function(ME, me.list) {
    pl <- list(spatial_me = me.list$spatial_me)
    if (me.list$spatial_me) {
        pl$prior_nux_true_alpha <- pl$prior_nux_true_beta <- vec.n.zeros(me.list$dx_me)
        if (inherits(ME$prior$car_rho, "NULL")) {
            lims <- 1 / range(ME$car_parts$lambda)
            pl$prior_rhox_true <- c(lims[1], lims[2])
        } else {
            pl$prior_rhox_true <- c(ME$prior$car_rho$lower, ME$prior$car_rho$upper)
        }        
    } else {
        pl$prior_rhox_true <- c(-1, 1)
        if (inherits(ME$prior$df, "NULL")) { 
            pl$prior_nux_true_alpha <- rep(3, times = me.list$dx_me)
            pl$prior_nux_true_beta <- rep(0.2, times = me.list$dx_me)
        } else {
            pl$prior_nux_true_alpha <- ME$prior$df$alpha
            pl$prior_nux_true_beta <- ME$prior$df$beta
        }
    }
    if (inherits(ME$prior$location, "NULL")) {
        pl$prior_mux_true_location   <- rep(0, times = me.list$dx_me)
        pl$prior_mux_true_scale      <- rep(100, times = me.list$dx_me)
    } else {
        pl$prior_mux_true_location   <- ME$prior$location$location
        pl$prior_mux_true_scale      <- ME$prior$location$scale
    }
    if (inherits(ME$prior$scale, "NULL")) {
        pl$prior_sigmax_true_df       <- rep(10, times = me.list$dx_me)
        pl$prior_sigmax_true_location <- rep(0, times = me.list$dx_me)
        pl$prior_sigmax_true_scale    <- rep(40, times = me.list$dx_me)        
    } else {
        pl$prior_sigmax_true_df       <- ME$prior$scale$df
        pl$prior_sigmax_true_location <- ME$prior$scale$location
        pl$prior_sigmax_true_scale    <- ME$prior$scale$scale
    }
    pl <- lapply(pl, as.array)
    pl$ME_prior_car_rho <- uniform(pl$prior_rhox_true[1], pl$prior_rhox_true[2])
    pl$ME_prior_df <- gamma(alpha = pl$prior_nux_true_alpha,
                            beta = pl$prior_nux_true_beta,
                            variable = names(ME$se)
                            )      
    pl$ME_prior_location <- normal(location = pl$prior_mux_true_location,
                               scale = pl$prior_mux_true_scale,
                               variable = names(ME$se)
                               )
    pl$ME_prior_scale <- student_t(df =  pl$prior_sigmax_true_df,
                                   location = pl$prior_sigmax_true_location,
                                   scale = pl$prior_sigmax_true_scale,
                                   variable = names(ME$se)
                                   )
    if (me.list$has_me) print_me_priors(pl, ME)
    return(pl)        
}

print_me_priors <- function(pl, ME) {
    if (pl$spatial_me) {
        if (inherits(ME$prior$car_rho, "NULL")) {
            message(
            "\n*Setting ME model prior",
            "\n*Setting prior for CAR spatial autocorrelation parameter (rho)"
            )
            print(pl$ME_prior_car_rho)
        }
    } else {
        if (inherits(ME$prior$df, "NULL")) {
            message(
            "\n*Setting ME model prior",
            "\n*Degrees of freedom (nu)"
            )
            print(pl$ME_prior_df)
        }
    }    
    if (inherits(ME$prior$location, "NULL")) {
        message(
            "\n*Setting ME model prior",
            "\n*Location (mean)"
            )
        print(pl$ME_prior_location)
    }
    if (inherits(ME$prior$scale, "NULL")) {
        message(
            "\n*Setting ME model prior",
            "\n*Scale"
        )
        print(pl$ME_prior_scale)
    }
}

