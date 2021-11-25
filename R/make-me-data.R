#' make_me_data
#' 
#' @description internal function to prepare Stan data list for observational error models
#' 
##' @param ME user-provided list of observational error stuff
##' 
##' @param x design matrix, excluding any Wx (slx) terms
##' 
##' @return Partial data list to pass to Stan, including observations, standard errors, and prior parameters.
##' 
##' @noRd
make_me_data <- function(ME, x) { # for x pass in x_no_Wx !
    x.df <- as.data.frame(x)
    n <- nrow(x.df)  
    if (inherits(ME, "NULL")) {
        me.list <- empty_me.list(x)
        return(me.list)
    }
    check_me_data(ME)
    if (!all(names(ME$se) %in% names(x.df))) stop("All column names in ME$se must be found in the model matrix (from model.matrix(formula, data)).")        
    # gather any/all variables without ME
    x_obs_idx <- as.array( which( !names(x.df) %in% names(ME$se) )) 
    x_obs <- as.data.frame(x.df[, x_obs_idx])
    dx_obs <- ncol(x_obs)
    # gather ME variables
    x_me_idx <- as.array(na.omit(match(names(ME$se), names(x.df))))
    x_me <- as.data.frame(x.df[,x_me_idx])
    dx_me <- ncol(x_me)
    sigma_me <- ME$se
    # check user-provided boundary condition
    if (any(x_me < ME$bounds[1]) || any(x_me > ME$bounds[2])) stop("In ME: bounded variable has elements outside of user-provided bounds (", ME$bounds[1], ", ", ME$bounds[2], ")")
    # handle condition when all covariates are modeled, none classified as perfectly 'observed'
    if (!dx_obs) {
        x_obs <- model.matrix(~ 0, x.df) 
        x_obs_idx <- a.zero()
    }
    me.list <- list(
        dx_obs = dx_obs,
        dx_me = dx_me,
        x_obs_idx = x_obs_idx,
        x_me_idx = x_me_idx,
        bounds = ME$bounds,
        x_obs = x_obs, 
        x_me = array(t(x_me), dim = c(dx_me, n)),
        sigma_me = array(t(sigma_me), dim = c(dx_me, n)),
        nm_me = names(ME$se),
        has_me = TRUE,
        spatial_me = ME$spatial_me,
        use_logit = as.array(as.numeric(ME$logit))
        )
    me.list <- c(me.list, ME$car_parts)
    pl <- make_me_priors(ME, me.list)
    me.list <- c(me.list, pl)
    return(me.list)
}

#' return items in data list ready for Stan when there is no ME model at all.
#' 
#' @param x design matrix, excluding any Wx (slx) terms
#' @noRd
empty_me.list <- function(x) {
    x.df <- as.data.frame(x)
    n <- nrow(x.df)  
    dx_me <- 0
    x_me_idx = a.zero()
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
        bounds = c(-Inf, Inf),          
        x_obs = x_obs,
        x_me = a.n.zeros(n),
        sigma_me = a.n.zeros(n),
        ME_prior_car_rho = vec.n.zeros(2),
        has_me = FALSE,
        spatial_me = FALSE,
        use_logit = as.array(rep(0, times = dx_me))
    )
    me.list <- c(me.list, car_parts_shell(n))
    pl <- make_me_priors(NULL, me.list)
    me.list <- c(me.list, pl)
    return(me.list)
} 

make_me_priors <- function(ME, me.list) {
    if (inherits(ME, "NULL")) {
        empty_priors <- vector(mode = "list", length = 8)
        names(empty_priors) <- c("prior_mux_true_location",
                                "prior_mux_true_scale",
                                "prior_sigmax_true_df",
                                "prior_sigmax_true_location",
                                "prior_sigmax_true_scale",
                                "prior_nux_true_alpha",
                                "prior_nux_true_beta",
                                "prior_rhox_true")
        for (i in 1:7) empty_priors[[i]] <- numeric(0)
        empty_priors[[8]] <- c(-1,1)
        empty_priors$spatial_me <- FALSE
        return (empty_priors)
    }
    pl <- list(spatial_me = me.list$spatial_me)    
    pl$prior_mux_true_location   <- as.array(ME$prior$location$location)
    pl$prior_mux_true_scale      <- as.array(ME$prior$location$scale)
    pl$prior_sigmax_true_df       <- as.array(ME$prior$scale$df)
    pl$prior_sigmax_true_location <- as.array(ME$prior$scale$location)
    pl$prior_sigmax_true_scale    <- as.array(ME$prior$scale$scale)
    if (me.list$spatial_me) {
        pl$prior_nux_true_alpha <- pl$prior_nux_true_beta <- a.zero()
        pl$prior_rhox_true <- c(ME$prior$car_rho$lower, ME$prior$car_rho$upper)
    } else {
        pl$prior_nux_true_alpha <- as.array(ME$prior$df$alpha)
        pl$prior_nux_true_beta <- as.array(ME$prior$df$beta)
        pl$prior_rhox_true <- c(-1, 1)        
    }
    return (pl)
}

