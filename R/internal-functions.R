#' Remove intercept from formula
#' 
#' @noRd
#' @param x A model frame or matrix with or without an intercept
#' @return The same model frame or matrix with any intercept term removed
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
}

#' scale model matrix
#'
#' @noRd
#' @param x model matrix
#' @return scaled model matrix. Binary indicator variables will not be scaled.
scale_x <- function(x, center, scale) {
    scale_params <- list()
    if (center) scale_params$center <- rep(0, times = ncol(x))
    if (scale) scale_params$scale <- rep(1, times = ncol(x))
    for (i in 1:ncol(x)) {
        l <- length(unique(x[,i]))
        if (l == 2) next
        x.tmp <- scale(x[,i], center = center, scale = scale)
        if (center) scale_params$center[i] <- attributes(x.tmp)$`scaled:center`
        if (scale) scale_params$scale[i] <- attributes(x.tmp)$`scaled:scale`
        x[,i] <- as.numeric(x.tmp)
    }
    return(list(x = x, params = scale_params))
    }

#' @importFrom stats model.matrix
#' @noRd
SLX <- function(f, DF, SWM) {
    # row-standardize the connectivity matrix
    W <- SWM / rowSums(SWM)
    x <- remove_intercept(model.matrix(f, DF))
    Wx <- W %*% x
    dimnames(Wx)[[2]] <- paste0("w.", dimnames(Wx)[[2]])
    return(Wx)
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
#' @param stanfit a stanfit object
#' @param par name of parameter
post_summary <- function(stanfit, par) {
  par_mat <- as.matrix(stanfit, par = par)
  par_summary <- data.frame(
    mean = apply(par_mat, 2, mean),
    sd = apply(par_mat, 2, sd),
    q.025 = apply(par_mat, 2, quantile, probs = .025),
    q.20 = apply(par_mat, 2, quantile, probs = .20),
    q.50 = apply(par_mat, 2, quantile, probs = .5),
    q.80 = apply(par_mat, 2, quantile, probs = .8),
    q.975 = apply(par_mat, 2, quantile, probs = .975)
  )
    return(par_summary)
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

#' Root mean square error
#' 
#' @noRd
#' @param error numeric vector of residuals
rmse <- function(error, digits) {
  rmse <- sqrt(mean(error^2))
  if(!missing(digits)) rmse <- round(rmse, digits = digits)
  return(rmse)
}

#' logit
#' @noRd
#' @param p probability, ratio
logit <- function(p) log(p/(1-p))

#' Build list of priors
#' @importFrom stats sd
#' @noRd
make_priors <- function(user_priors = NULL, y, x, xcentered, rhs_scale_global, scaling_factor = 2, link = c("identity", "log", "logit"), EV, offset) {
  if (link == "identity") scale.y <- sd(y) else scale.y <- 1
  if (link == "log") {
      if (any(y == 0)) y[which(y == 0)] <- 1 # local assignment only, not returned
      y <- log(y / exp(offset))
      scale.y <- sd(y)
      }
  alpha_scale <- scaling_factor * scale.y
  if (xcentered) alpha_mean <- mean(y) else alpha_mean <- 0
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
    beta_location <- rep(0, times = ncol(x))
    priors$beta <- cbind(beta_location, beta_scale)
    dimnames(priors$beta)[[1]] <- dimnames(x)[[2]]
  } else {
      priors$beta <- matrix(0, nrow = 0, ncol = 2)
  }
  # the following will be ignored when when not needed (RE scale, resid scale, student T df)
  priors$alpha_tau <- c(df = 20, location = 0, scale = scaling_factor * scale.y)
  priors$sigma <- c(df = 10, location = 0, scale = scaling_factor * scale.y)
  priors$nu <- c(alpha = 3, beta = 0.2)
  if ("rhs" %in% names(user_priors)) {
      scale_ev <- sd(EV[,1])
      priors$rhs = c(slab_df = 15, slab_scale = 0.5 * (scale.y / scale_ev), scale_global = rhs_scale_global)
      }
  for (i in seq_along(user_priors)) {
      if (!is.null(user_priors[[i]])) {
          par <- names(user_priors)[i]
          priors[[which(names(priors) == par)]] <- user_priors[[i]]
      }      
  }
  return(priors)
}

#' Log sum of exponentials
#' @noRd
log_sum_exp <- function(x) {
  # learned/stolen from R. McElreath's Rethinking package
  xmax <- max(x)
  xsum <- sum( exp( x - xmax ) )
  xmax + log(xsum)
}

IDW <- function(mat, lambda) {
    mat <- 1/mat^lambda
    diag(mat) <- 0
    mat <- mat / rowSums(mat)
    class(mat) <- "matrix"
    return(mat)
}

clean_results <- function(samples, pars, is_student, has_re, C, Wx, x, x_me_unbounded_idx, x_me_bounded_idx) {
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
    if (sum(x_me_unbounded_idx)) {
        x_names <- dimnames(x)[[2]]        
        for (i in seq_along(x_me_unbounded_idx)) {
            x.id <- paste0("x_", rep(x_names[x_me_unbounded_idx[i]], times = n), paste0("[", 1:n, "]"))
            names(samples)[grep(paste0("x_true_unbounded\\[", i, ","), names(samples))] <- x.id          
        }
    }   
    if (sum(x_me_bounded_idx)) {
        x_names <- dimnames(x)[[2]]
        for (i in seq_along(x_me_bounded_idx)) {        
            x.id <- paste0("x_", rep(x_names[x_me_bounded_idx[i]], times = n), paste0("[", 1:n, "]"))
            names(samples)[grep(paste0("x_true_bounded\\[", i, ","), names(samples))] <- x.id
        }
    }
    if ("sigma" %in% pars) samples <- par_alias(samples, "^sigma\\[1\\]", "sigma")
    if ("rho" %in% pars) samples <- par_alias(samples, "^rho\\[1\\]", "rho")
    if ("theta_scale" %in% pars) samples <- par_alias(samples, "^theta_scale\\[1\\]", "theta_scale")       
    if (is_student) samples <- par_alias(samples, "^nu\\[1\\]", "nu")
    if (has_re) samples <- par_alias(samples, "^alpha_tau\\[1\\]", "alpha_tau")
    main_pars <- pars[which(pars %in% c("nu", "intercept", "alpha_tau", "gamma", "beta", "sigma", "rho", "spatial_scale", "theta_scale"))]
    summary_list <- lapply(main_pars, post_summary, stanfit = samples)
    summary <- do.call("rbind", summary_list)
    summary <- apply(summary, 2, round, 3)
  residuals <- as.matrix(samples, pars = "residual")
    residuals <- apply(residuals, 2, mean)
  if (is.logical(C)) {
      Residual_MC <- C <- Expected_MC <- NA
  } else {
      Residual_MC <- round(mc(residuals, w = C), 3)
      if (!has_b) {
          n <- length(residuals)
          Expected_MC <- -1/(n - 1)
      } else {
          Expected_MC <- expected_mc(X = x, C = C)
      }
     }
    WAIC <- geostan::waic(samples)
    diagnostic <- c(WAIC = as.numeric(WAIC[1]), Eff_pars = as.numeric(WAIC[2]), Lpd = as.numeric(WAIC[3]),
                    Residual_MC = Residual_MC, Expected_MC = Expected_MC)
    out <- list(summary = summary, diagnostic = diagnostic, stanfit = samples)
    return(out)
}

#' Print priors that were set automatically
#' 
#' @noRd
#' @param user_priors list of prior specifications provided by user, with missing priors set to NULL
#' @param priors list of priors returned by make_priors, after dropping unneeded parameters.
#' @return Prints out information using \code{message()}
print_priors <- function(user_priors, priors) {
    printed = 0
  for (i in seq_along(user_priors)) {
      nm <- names(user_priors)[i]
      u.p. <- user_priors[[i]]
      if (is.null(u.p.) & nm %in% names(priors)) {          
          message("\n*Setting prior parameters for ", nm)
          p <- priors[[which(names(priors) == nm)]]
          if (nm == "beta") {
              message("Gaussian")
              colnames(p) <- c("Location", "Scale")
              print(p, row.names = FALSE)
          }          
          if (nm == "intercept") {
              message("Gaussian")
              message("Location: ", p[1])
              message("Scale: ", p[2])
          }
          if (nm == "sigma") {
              message("Student's t")
              message("Degrees of freedom: ", p[1])
              message("Location: ", p[2])
              message("Scale: ", p[3])
          }
          if (nm == "nu") {
              message("Gamma")
              message("alpha: ", p[1])
              message("beta: ", p[2])
          }
          if (nm == "alpha_tau") {
              message("Student's t")
              message("Degrees of freedom: ", p[1])
              message("Location: ", p[2])
              message("Scale: ", p[3])
          }
          if (nm == "rhs") {
              message("Global shrinkage prior (scale_global): ", p["scale_global"])
              message("Slab degrees of freedom: ", p["slab_df"])
              message("Slab scale: ", p["slab_scale"])
          }
  }  
  } 
}

