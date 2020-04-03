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

#' scale model matrix
#'
#' @noRd
#' @param x model matrix
#' @return scaled model matrix. Binary indicator variables will not be scaled.
scale_x <- function(x, center, scale) {
    for (i in seq_along(ncol(x))) {
        l <- length(unique(x[,i]))
        if (l == 2) next
        x[,i] <- as.numeric(scale(x[,i], center = center, scale = scale))
    }
    return(x)
    }

#' @importFrom stats model.matrix
#' @noRd
SLX <- function(f, DF, SWM, cx, sx) {
    # row-standardize the connectivity matrix
    W <- SWM / rowSums(SWM)
    x <- remove_intercept(model.matrix(f, DF))
    x <- apply(x, MARGIN = 2, FUN = scale, center = cx, scale = sx)
    Wx <- W %*% x
    dimnames(Wx)[[2]] <- paste0("sl.", dimnames(Wx)[[2]])
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
    o <- model.offset(model.frame(f, df))
    if (is.null(o)) return(cbind(y, x))
    return(cbind(y, o, x))
}

#' Coerce a vector to data frame with the orignal values and integer ID values
#' 
#' @noRd
#' @param id a vector
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

#' 
#' @noRd
family_2_integer <- function(family) {
    if (family == "gaussian") return(1)
    if (family == "student_t") return(2)
    if (family == "binomial") return(3)
    if (family == "poisson") return(4)
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
make_priors <- function(user_priors = NULL, y, x, xcentered, rhs_scale_global, scaling_factor = 2.5, link = c("identity", "log", "logit"), EV) {
  if (link == "log") y <- log(y)
  if (link == "logit") y <- logit(y[,1] / (y[,1] + y[,2]))
  scaley <- sd(y)
  alpha_scale <- max(10 * sd(y), 5)
  alpha_mean <- 0
  alpha <- c(df = 15, location = alpha_mean, scale = alpha_scale)
  alpha <- round(alpha)
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
    beta_scale <- scaling_factor * (scaley / scalex)
    beta_location <- rep(0, times = ncol(x))
    beta_df <- rep(15, times = ncol(x))
    priors$beta <- cbind(beta_df, beta_location, beta_scale)
    priors$beta <- round(priors$beta, 4)
    dimnames(priors$beta)[[1]] <- dimnames(x)[[2]]
  } else {
      priors$beta <- matrix(0, nrow = 0, ncol = 3)
  }
  # the following will be ignored when when not needed (RE scale, resid scale, student T df)
  priors$alpha_tau <- round(c(df = 20, location = 0, scale = scaling_factor * scaley))
  priors$sigma <- round(c(df = 5, location = 0, scale = scaling_factor * scaley))
  priors$nu <- c(alpha = 2, beta = .1)
  if ("rhs" %in% names(user_priors)) {
      scale.ev <- sd(EV[,1])
      priors$rhs = c(slab_df = 15, slab_scale = scaling_factor * (scaley / scale.ev), scale_global = rhs_scale_global)
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

clean_results <- function(samples, pars, is_student, has_re, C, x) {
  has_b <- "beta" %in% pars  
  if (has_b) {
    b_names = paste0("b_", dimnames(x)[[2]])
    samples <- par_alias(samples, "^beta\\[", b_names)
  }
  if ("sigma[1]" %in% pars) samples <- par_alias(samples, "^sigma\\[1\\]", "sigma")
  if (is_student) samples <- par_alias(samples, "^nu\\[1\\]", "nu")
  if (has_re) samples <- par_alias(samples, "^alpha_tau\\[1\\]", "alpha_tau")
  main_pars <- pars[which(pars %in% c("intercept", "alpha_tau", "beta", "sigma", "nu", "rho"))]
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
  RMSE <- rmse(residuals, digits = 2)
  WAIC <- geostan::waic(samples)
  diagnostic <- c(WAIC = as.numeric(WAIC[1]), Eff_pars = as.numeric(WAIC[2]), Lpd = as.numeric(WAIC[3]),
                  Residual_MC = Residual_MC, Expected_MC = Expected_MC, RMSE = RMSE)
  out <- list(summary = summary, diagnostic = diagnostic, stanfit = samples)
  return(out)
}

