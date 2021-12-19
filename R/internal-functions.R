#' Remove intercept from formula
#' 
#' @noRd
#' @param x A model frame or matrix with or without an intercept
#' @return The same model frame or matrix but without the intercept term 
remove_intercept <- function(x) {
  varnames <- dimnames(x)[[2]][-which(attributes(x)$assign == 0)]
  xnew <- as.matrix(x[,which(attributes(x)$assign != 0)])
  dimnames(xnew)[[2]] <- varnames
  return(xnew)
}

#' Convert family to integer
#'
#' @noRd
#' @param family class family
#' @return interger value
#' 
family_2_int <- function(family) {
    if (family$family == "gaussian") return (1)
    if (family$family == "student_t") return (2)
    if (family$family == "poisson") return (3)
    if (family$family == "binomial") return (4)
    if (family$family == "auto_gaussian") return (5)
}


#' @noRd
#' @param standata The list of data as passed to rstan::sampling
#' @param samples The stanfit object returned by rstan::sampling
#' If centerx = TRUE, return the values on which covariates were centered. Handle ME variables appropriately using modeled covariate mean.
get_x_center <- function(standata, samples) {
    if (standata$center_x == 0) return (FALSE)
    dx_obs <- standata$dx_obs
    dx_me <- standata$dx_me
    M <- dx_obs + dx_me
    x_center <- vector(mode = "numeric", length = M)    
    if (dx_obs > 0) {
        x_obs_mean <- Matrix::colMeans(standata$x_obs)
        x_center[standata$x_obs_idx] <- x_obs_mean
    }
    if (dx_me > 0) {
        x_true <- as.matrix(samples, pars = 'mu_x_true')
        x_true_mean <- apply(x_true, 2, mean)
        ulo <- standata$use_logit
        for (k in 1:dx_me) if (ulo[k]) x_true_mean[k] <- inv_logit(x_true_mean[k])     
        x_center[standata$x_me_idx] <- x_true_mean
    }
    return (x_center)
}
    
#' Create matrix of spatially lagged covariates, using prepared (i.e., possibly centered) model matrix.
#' 
#' @importFrom stats model.matrix
#'
#' @param f formula (from slx)
#' @param DF `as.data.frame(data)`
#' @param x The possibly centered model matrix, no intercept
#' @param W row-standardized spatial weights matrix
#' 
#' @noRd
#' 
SLX <- function(f, DF, x, W) {
    z <- remove_intercept(model.matrix(f, DF))    
    x_slx <- as.matrix(x[, colnames(z)])
    Wx <- W %*% x_slx
    colnames(Wx) <- paste0("w.", colnames(z))
    return (Wx)       
   }

#' Prep data to return to user
#'
#' @noRd
#' @param frame Model frame from calling model.frame(formula, data, na.action = NULL)
#' @param x Model matrix, including SLX terms if any, after centering and scaling
#' @param y_mis_idx Index of censored counts in the outcome, if any.
make_data <- function(frame, x, y_mis_idx) {
    y <- model.response(frame)
    if (length(y_mis_idx) > 0) y[y_mis_idx] <- NA
    offset <- model.offset(frame)
    if (is.null(offset)) return(cbind(y, x))
    return(cbind(y, offset, x))
}

#' Coerce a vector to data frame with the orignal values and integer ID values
#' 
#' @noRd
#' @param id a vector
#' @param n scalar value, number of zero observations to fill.
#' @return a data.frame with columns id and idx
to_index <- function(id, n) {
  if(length(id) == 1) {
    return(data.frame(id = rep(0, times = n), idx = 0))
  }
  n_ids <- length(unique(id))
  idx <- data.frame(id = unique(id), idx = 1:n_ids)
  res <- merge(data.frame(rowid = 1:length(id), id = id), idx, by = "id")
  res[order(res$rowid), c("id", "idx")]
}

#' Summarize samples from an geostan_fit object
#' @import stats
#' @noRd
#' @param samples matrix of MCMC samples
#' 
post_summary <- function(samples) {
    qs = apply(samples, 2, quantile, probs = c(0.025, 0.2, 0.50, 0.80, 0.975))
    ms <- apply(samples, 2, mean)
    sds <- apply(samples, 2, sd)
    x <- as.data.frame(cbind(mean = ms, sd = sds, t(qs)))
    return (x)
}

#' logit
#' @noRd
#' @param p probability
logit <- function(p) log(p/(1-p))

#' inverse logit
#' @noRd
inv_logit <- function(x) exp(x)/(1 + exp(x))

#' Build list of priors
#' 
#' @importFrom stats sd
#' @noRd
make_priors <- function(user_priors = NULL, y, x, hs_global_scale, scaling_factor = 2, link = c("identity", "log", "logit"), EV, offset) {
    link <- match.arg(link)
    if (link == "identity") {
        scale.y <- sd(y)
        alpha_scale <- max(4 * sd(y), 5)
        alpha_mean <- mean(y)
    }
  if (link == "log") {
      if (any(y == 0)) y[which(y == 0)] <- 1 # local assignment only, not returned
      y <- log(y / exp(offset))
      alpha_mean <- mean(y)
      scale.y <- sd(y)
      alpha_scale <- max(4 * scale.y, 5)
  }
  if (link == "logit") {
      y <- y[,1] / (y[,1] + y[,2])
      alpha_mean <- 0
      scale.y <- sd(y)
      alpha_scale <- 5
  }
    priors <- list()
    priors$intercept <- normal(location = alpha_mean, scale = alpha_scale)
  if (ncol(x)) {
    scalex <- vector(mode = "numeric", length = ncol(x))
    nxvals <- apply(x, 2, FUN = function(x) length(unique(x)))
    x_cont <- which(nxvals > 2)
    x_bin <- which(nxvals == 2)
    scalex_cont <- apply(as.matrix(x[,x_cont]), 2, sd)
    scalex[x_cont] <- scalex_cont
    if (length(x_bin)) {
      scalex_bin <- apply(as.matrix(x[,x_bin]), 2, function(x) max(x) - min(x))
      scalex[x_bin] <- scalex_bin
    }
    beta_scale <- scaling_factor * (scale.y / scalex)
    for (i in 1:length(beta_scale)) beta_scale[i] <- max(beta_scale[i], 5)
    beta_location <- rep(0, times = ncol(x))
    priors$beta <- normal(location = beta_location, scale = beta_scale, variable = colnames(x))
  } else {
      priors$beta <- normal(location = 0, scale = 1)
  }
  priors$alpha_tau <- student_t(df = 10, location = 0, scale = max(scaling_factor * scale.y, 3))
  priors$sigma <- student_t(df = 10, location = 0, scale = max(scaling_factor * scale.y, 3))
  priors$nu <- gamma(alpha = 3, beta = 0.2)
  if(!missing(hs_global_scale)) {
      scale_ev <- sd(EV[,1])
      priors$beta_ev = hs(global_scale = hs_global_scale, slab_df = 15, slab_scale = 0.5 * (scale.y / scale_ev))
  }
  for (i in seq_along(priors)) {
      par.name <- names(priors)[i]
      if(!is.null(user_priors[[par.name]])) {
          priors[[i]] <- user_priors[[par.name]]
      }      
  }
  return(priors)
}

#' Log sum of exponentials
#' @noRd
#'
#' @details Code adapted from Richard McElreath's Rethinking package, and other sources.
#' 
log_sum_exp <- function(x) {
  xmax <- max(x)
  xsum <- sum( exp( x - xmax ) )
  xmax + log(xsum)
}

clean_results <- function(samples, pars, is_student, has_re, Wx, x, x_me_idx) {
    n <- nrow(x)
    if ("gamma" %in% pars) {
    g_names = dimnames(Wx)[[2]]
    samples <- par_alias(samples, "^gamma\\[", g_names)        
    }
    has_b <- "beta" %in% pars  
    if (has_b) {  
        b_names = dimnames(x)[[2]] 
        samples <- par_alias(samples, "^beta\\[", b_names)
    }
    if (sum(x_me_idx)) {
        x_names <- dimnames(x)[[2]]        
        for (i in seq_along(x_me_idx)) {
            x.id <- paste0("x_", rep(x_names[x_me_idx[i]], times = n), paste0("[", 1:n, "]"))
            samples <- par_alias(samples, paste0("^x_true\\[" , i, ","), x.id)
        }
    }   
    if ("sigma" %in% pars) samples <- par_alias(samples, "^sigma\\[1\\]", "sigma")
    if ("rho" %in% pars) samples <- par_alias(samples, "^rho\\[1\\]", "rho")
    if ("theta_scale" %in% pars) samples <- par_alias(samples, "^theta_scale\\[1\\]", "theta_scale")
    if ("spatial_scale" %in% pars) samples <- par_alias(samples, "^spatial_scale\\[1\\]", "spatial_scale")
    if ("car_rho" %in% pars) {
        samples <- par_alias(samples, "^car_rho\\[1\\]", "car_rho")
        samples <- par_alias(samples, "^car_scale\\[1\\]", "car_scale")
    }
    if (is_student) samples <- par_alias(samples, "^nu\\[1\\]", "nu")
    if (has_re) samples <- par_alias(samples, "^alpha_tau\\[1\\]", "alpha_tau")
    main_pars <- pars[which(pars %in% c("nu", "intercept", "alpha_tau", "gamma", "beta", "sigma", "rho", "spatial_scale", "theta_scale", "car_scale", "car_rho"))]
    S <- as.matrix(samples, pars = main_pars)
    summary <- post_summary(S)
    Residual_MC <- NA
    if ("log_lik" %in% pars) WAIC <- geostan::waic(samples) else WAIC <- rep(NA, 3)
    diagnostic <- c(WAIC = as.numeric(WAIC[1]), Eff_pars = as.numeric(WAIC[2]), Lpd = as.numeric(WAIC[3]),
                    Residual_MC = Residual_MC)
    out <- list(summary = summary,
                diagnostic = diagnostic, stanfit = samples)
    class(out) <- append("geostan_fit", class(out))
    return(out)
}

#' Rename Stan model parameters
#' 
#' @noRd
#' @param samples stan model object
#' @param par original parameter name to search for; regular expression
#' @param replacement new paramter name or names
par_alias <- function(samples, par, replacement) {
  found <- grep(par, names(samples))
  if(length(found) != length(replacement)) message("Attemping to rename parameter ", par, "; replacement is of wrong length. Is of length ", length(replacement), ", must be ", length(found))
  chains <- length(samples@sim$samples)
  for (i in 1:chains) names(samples@sim$samples[[i]])[ grep(par, names(samples@sim$samples[[i]])) ] <- replacement
  names(samples)[ grep(par, names(samples)) ] <- replacement
  return(samples)
}

#' Print priors that were set automatically
#' 
#' @noRd
#' @param user_priors list of prior specifications provided by user, with missing priors set to NULL
#' @param priors list of priors returned by make_priors, after dropping unneeded parameters.
#' @return Prints out information using \code{message()}
print_priors <- function(user_priors, priors) {
  for (i in seq_along(priors)) {
      nm <- names(priors)[i]
      u.p. <- user_priors[[nm]]
      if (is.null(u.p.)) {          
          p <- priors[[nm]]
          if (nm == "intercept") {
              message("\n*Setting prior parameters for ", nm)              
              #message("Gaussian")
              print(p)
          }          
          if (nm == "beta") {
              message("\n*Setting prior parameters for ", nm)              
              print(p)
          }          
          if (nm == "sigma") {
              message("\n*Setting prior parameters for ", nm)              
              print(p)
          }
          if (nm == "nu") {
              message("\n*Setting prior parameters for ", nm)              
              print(p)
          }
          if (nm == "alpha_tau") {
              message("\n*Setting prior parameters for ", nm)              
              print(p)
          }
          if (nm == "beta_ev") {
              message("\n*Setting prior for eigenvector coefficients")              
              print(p)
          }
          if (nm == "car_scale") {
              message("\n*Setting prior for CAR scale parameter (car_scale)")
              print(p)
          }
          if (nm == "car_rho") {
              message("\n*Setting prior for CAR spatial autocorrelation parameter (rho)")
              print(p)
          }
      }  
  } 
}


a.zero <- function() array(0, dim = 1)
a.n.zeros <- function(n) array(0, dim = c(0, n))
vec.n.zeros <- function(n) rep(0, n)

#' return empty car_parts list
#' @noRd
#' @param n vector length
car_parts_shell <- function(n) {
    list(
      nC = 1,
      nAx_w = 1,
      C = array(1, dim = c(1, 1)),
      Delta_inv = vec.n.zeros(n),
      log_det_Delta_inv = 0,
      Ax_w = a.zero(),
      Ax_v = a.zero(),
      Ax_u = vec.n.zeros(n+1), 
      Cidx = a.zero(),
      lambda = vec.n.zeros(n),
      WCAR = 0
    )
}

#' @noRd
empty_icar_data <- function(n) {
    dl <- list(
        type = 0,
        k = 1,
        m = 0,
        group_size = a.zero(),
        group_idx = vec.n.zeros(n),
        A = model.matrix(~ 0, data.frame(a=1:n)),
        n_edges = 1,
        node1= as.array(1),
        node2= as.array(1),
        weight = as.array(1),
        comp_id = as.array(rep(1, n)),
        inv_sqrt_scale_factor = as.array(1)
    )
    return (dl)
}

#' @noRd
empty_esf_data <- function(n) {
    list(
        dev = 0,
        EV = model.matrix(~ 0 , data.frame(a=1:n)),
        global_scale = 0,
        slab_scale = 0,
        slab_df = 0        
    )
}

#' @noRd
empty_car_data <- function() {
    list(
        car = 0,
        car_rho_lims = c(-1, 1)
    )
}
