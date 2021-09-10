

#' check data (x, w) dimensions and class for spatial autocorrelation measures
#' @noRd
check_sa_data <- function(x, w) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix") | inherits(w, "Matrix"))
    stopifnot(all(dim(w) == length(x)))
}


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

#' center model matrix
#'
#' @noRd
#' 
#' @param x model matrix
#' 
#' @param center Logical indicator for centering covariates, or numeric vector of values on which covariates will be centered.s
#' 
#' @return centered model matrix. No special hanlding of indicator variables.
#'
center_x <- function(x, center) {
    x <- scale(x, center = center, scale = FALSE)
    x_center <- attributes(x)$`scaled:center`
    if (center[1] || length(center) > 1) {
        message("Centering covariates:")
        for (i in 1:ncol(x)) message(" ", colnames(x)[i], ": ",  x_center[i])
    }
    return(x)
}

#' Create matrix of spatially lagged covariates, using prepared (i.e., possibly centered) model matrix.
#' 
#' @importFrom stats model.matrix
#'
#' @param f formula (from slx)
#'
#' @param DF `as.data.frame(data)`
#'
#' @param x The possibly centered model matrix, no intercept
#'
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
#' @param f formula
#' @param x Model matrix, including SLX terms if any, after centering and scaling
#' @param df user provided data.frame
make_data <- function(f, df, x) {
    y <- model.response(model.frame(f, df))
    offset <- model.offset(model.frame(f, df))
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

#' Rename Stan model parameters
#' 
#' @noRd
#' @param samples stan model object
#' @param par original parameter name to search for; regular expression
#' @param replacement new paramter name or names
par_alias <- function(samples, par, replacement) {
  found <- grep(par, names(samples))
  if(length(found) != length(replacement)) message("replacement is of wrong length. Is of length ", length(replacement), ", must be ", length(found),
                                                "\nAttempt to rename parameter ", par, " was allowed but may cause problems.")
  names(samples@sim$samples[[1]])[ grep(par, names(samples@sim$samples[[1]])) ] <- replacement
  names(samples)[ grep(par, names(samples)) ] <- replacement
  return(samples)
}

#' logit
#' @noRd
#' @param p probability
logit <- function(p) log(p/(1-p))

#' Build list of priors
#' @importFrom stats sd
#' @noRd
make_priors <- function(user_priors = NULL, y, x, rhs_scale_global, scaling_factor = 2, link = c("identity", "log", "logit"), EV, offset) {
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
  alpha <- c(location = alpha_mean, scale = alpha_scale)
  priors <- list(intercept = alpha)
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
    priors$beta <- cbind(beta_location, beta_scale)
    dimnames(priors$beta)[[1]] <- dimnames(x)[[2]]
    stopifnot(ncol(x) == nrow(priors$beta))      
  } else {
      priors$beta <- matrix(0, nrow = 0, ncol = 2)
  }
  # the following will be ignored when when not needed (RE scale, resid scale, student T df)
  priors$alpha_tau <- c(df = 10, location = 0, scale = max(scaling_factor * scale.y, 3))
  priors$sigma <- c(df = 10, location = 0, scale = max(scaling_factor * scale.y, 3))
  priors$nu <- c(alpha = 3, beta = 0.2)
  if(!missing(rhs_scale_global)) {
      scale_ev <- sd(EV[,1])
      priors$rhs = c(slab_df = 15, slab_scale = 0.5 * (scale.y / scale_ev), scale_global = rhs_scale_global)
      }
  for (i in seq_along(priors)) {
      par.name <- names(priors)[i]
      if(!is.null(user_priors[[par.name]])) {
          if (par.name == "beta") {
              beta_prior <- user_priors[[par.name]]
              if (ncol(x) == 1 &
                  inherits(beta_prior, "numeric")
                  ) user_priors[[par.name]] <- beta_prior <- as.data.frame(matrix(beta_prior, ncol = 2))                  
              stopifnot(ncol(x) == nrow(beta_prior))      
          }          
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
            names(samples)[grep(paste0("x_true\\[", i, ","), names(samples))] <- x.id          
        }
    }   
    if ("sigma" %in% pars) samples <- par_alias(samples, "^sigma\\[1\\]", "sigma")
    if ("rho" %in% pars) samples <- par_alias(samples, "^rho\\[1\\]", "rho")
    if ("theta_scale" %in% pars) samples <- par_alias(samples, "^theta_scale\\[1\\]", "theta_scale")
    if ("spatial_scale" %in% pars) samples <- par_alias(samples, "^spatial_scale\\[1\\]", "spatial_scale")
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
          if (nm == "beta") {
              message("\n*Setting prior parameters for ", nm)              
              message("Gaussian")
              colnames(p) <- c("Location", "Scale")
              print(p, row.names = FALSE)
          }          
          if (nm == "intercept") {
              message("\n*Setting prior parameters for ", nm)              
              message("Gaussian")
              print(p)
          }
          if (nm == "sigma") {
              message("\n*Setting prior parameters for ", nm)              
              message("Student's t")
              print(p)
          }
          if (nm == "nu") {
              message("\n*Setting prior parameters for ", nm)              
              message("Gamma")
              message("alpha: ", p[1])
              message("beta: ", p[2])
          }
          if (nm == "alpha_tau") {
              message("\n*Setting prior parameters for ", nm)              
              message("Student's t")
              print(p)
          }
          if (nm == "rhs") {
              message("\n*Setting prior parameters for eigenvector coefficients")              
              message("Regularized horseshoe (RHS) prior parameters:")
              message("Global shrinkage prior (scale_global): ", p["scale_global"])
              message("Slab degrees of freedom: ", p["slab_df"])
              message("Slab scale: ", p["slab_scale"])
          }
          if (nm == "car_scale") {
              message("\n*Setting prior parameters for car_scale")
              message("Student's t")
              print(p)

          }
  }  
  } 
}



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

a.zero <- function() array(0, dim = 1)
a.n.zeros <- function(n) array(0, dim = c(0, n))
vec.n.zeros <- function(n) rep(0, n)

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
        scale_global = 0,
        slab_scale = 0,
        slab_df = 0
    )
}
