#' Spatial autocorrelation estimator
#'
#' @description The approximate-profile likelihood estimator for the spatial autocorrelation parameter from a simultaneous autoregressive (SAR) model.
#' @param x Numeric vector of values, length n. This will be standardized internally with \code{scale(x)}.
#' @param w An n x n row-standardized spatial connectivity matrix. See \link[geostan]{shape2mat}.
#' @param digits Number of digits to round results to; defaults to \code{digits = 3}.
#' @return the APLE estimate.
#' @export
#' @source
#'
#' Li, Honfei and Calder, Catherine A. and Cressie, Noel (2007). Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis: 39(4): 357-375.
#'
#' @seealso \link[geostan]{mc} \link[geostan]{moran_plot} \link[geostan]{lisa} \link[geostan]{sim_sar}
#' @examples
#' 
#' library(sf)
#' data(ohio)
#' w <- shape2mat(ohio, "W")
#' x <- ohio$unemployment
#' aple(x, w)
#'
aple <- function(x, w, digits = 3) {
    z <- as.numeric(scale(x))
    n <- length(z)
    I <- diag(n)
    lambda <- eigen(w)$values
    w2 <- (w + t(w)) / 2
    top <- t(z) %*% w2 %*% z
    wl <- (t(w) %*% w + as.numeric(t(lambda) %*% lambda) * I / n)
    bottom <- t(z) %*% wl %*% z
    round(as.numeric( top / bottom ), digits = digits)
}

#' Simulate spatially autocorrelated data
#'
#' @description Given a spatial weights matrix and degree of autocorrelation, returns autocorrelated data. 
#' @export
#' @param n The number of samples required. Defaults to \code{n=1} to return a \code{k}-length vector; if \code{n>1}, an \code{n x k} matrix is returned (i.e. each row will contain a sample of correlated values).
#' @param mu A \code{k}-length vector of mean values. Defaults to a vector of zeros with length equal to \code{nrow(W)}.
#' @param w Row-standardized \code{k x k} spatial weights matrix.
#' @param rho Spatial autocorrelation parameter in the range [-1, 1]. Typically a scalar value; otherwise a K-length numeric vector.
#' @param sigma Scale parameter (standard deviation). Defaults to \code{sigma = 1}. Typically a scalar value; otherwise a K-length numeric vector.
#' @param ... further arguments passed to \code{MASS::mvrnorm}.
#' 
#' @return
#'
#' If \code{n = 1} a vector of the same length as \code{mu}, otherwise an \code{n x length(mu)} matrix with one sample in each row.
#'
#' @details Calls \code{MASS::mvrnorm} internally to draw from the multivariate normal distribution. The covariance matrix is specified following the simultaneous autoregressive (SAR) model. 
#'
#' @seealso \link[geostan]{aple} 
#' 
#' @examples
#' 
#' data(ohio)
#' w <- shape2mat(ohio, "W")
#' x <- sim_sar(w=w, rho=.8)
#' aple(x, w)
#'
#' @seealso \link[geostan]{shape2mat} \link[MASS]{mvrnorm} 
#'
#' @importFrom MASS mvrnorm
#' 
sim_sar <- function(n = 1, mu = rep(0, nrow(w)), w, rho, sigma = 1, ...) {
    if (!inherits(w, "matrix") | mode(w) != "numeric" | nrow(w) != ncol(w) | !all(rowSums(w) %in% c(0,1))) stop("W must be a square, row-standardized numeric matrix.")
    K <- nrow(w)
    if (missing(mu)) {
        mu <- rep(0, K)
    } else {
        if (length(mu) != K | !inherits(mu, "numeric")) stop("mu must be a numeric vector with length equal to nrow(W).")
        }
    if (!inherits(rho, "numeric") | !length(rho) %in% c(1, K) | any(rho > 1 | rho < -1)) stop("rho must be numeric value within range [-1, 1], or a k-length numeric vector where K=nrow(W).")
    if (!inherits(sigma, "numeric") | !length(sigma) %in% c(1, K) | !all(sigma > 0)) stop("sigma must be a positive numeric value, or k-length numeric vector, with K=nrow(W).")
    I <- diag(K)
    S <- crossprod(solve(I - rho * w)) * sigma
    x <- MASS::mvrnorm(n = n, mu = mu, Sigma = S, ...)
    return(x)
}

#' The Moran coefficient
#'
#' @description The Moran coefficient, a measure of spatial autocorrelation (also known as Global Moran's I)
#' @export
#' @param x Numeric vector of input values, length n.
#' @param w An n x n spatial connectivity matrix. See \link[geostan]{shape2mat}. 
#' @param digits Number of digits to round results to; defaults to \code{digits = 3}.
#' @return The Moran coefficient, a numeric value.
#'
#' @details If any observations with no neighbors are found (i.e. \code{any(rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#'
#' @seealso \link[geostan]{moran_plot} \link[geostan]{lisa} \link[geostan]{aple}
#' 
#' @examples
#' 
#' library(sf)
#' data(ohio)
#' w <- shape2mat(ohio, style = "W")
#' x <- ohio$unemployment
#' mc(x, w)
#'
#'
mc <- function(x, w, digits = 3) {
    if(missing(x) | missing(w)) stop("Must provide data x (length n vector) and n x n spatial weights matrix (w).")    
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
        message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data.")
        x <- x[-zero.idx]
        w <- w[-zero.idx, -zero.idx]
    }
    xbar <- mean(x)
    z <- x - xbar
    ztilde <- as.numeric(w %*% z)
    A <- sum(rowSums(w))
    n <- length(x)
    mc <- as.numeric( n/A * (z %*% ztilde) / (z %*% z))
    return(round(mc, digits = digits))
}

#' Moran plot
#'
#' @description Plots a set of values against their spatially lagged values and gives the Moran coefficient as a measure of spatial autocorrelation. 
#' @export
#' @import ggplot2
#' @param y A numeric vector of length n.
#' @param w An n x n spatial connectivity matrix.
#' @param xlab Label for the x-axis. 
#' @param ylab Label for the y-axis.
#' @param pch Symbol type.
#' @param col Symbol color. 
#' @param size Symbol size.
#' @param alpha Symbol transparency.
#' @param lwd Width of the regression line. 
#' @details For details on the symbol parameters see the documentation for \link[ggplot2]{geom_point}.
#'
#' If any observations with no neighbors are found (i.e. \code{any(rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#' 
#' @return Returns a \code{gg} plot, a scatter plot with \code{y} on the x-axis and its spatially lagged values on the y-axis (i.e. a Moran plot).
#'
#' @seealso \link[geostan]{mc} \link[geostan]{lisa} \link[geostan]{aple}
#' 
#' @examples
#' 
#' library(sf)
#' data(ohio)
#' y <- ohio$unemployment
#' w <- shape2mat(ohio, "W")
#' moran_plot(y, w)
#'
moran_plot <- function(y, w, xlab = "y (centered)", ylab = "Spatial Lag", pch = 20, col = "darkred", size = 2, alpha = 1, lwd = 0.5) {
    if (!(inherits(y, "numeric") | inherits(y, "integer"))) stop("y must be a numeric or integer vector")
    sqr <- all( dim(w) == length(y) )
    if (!inherits(w, "matrix") | !sqr) stop("w must be an n x n matrix where n = length(y)")
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
        message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data.")
        y <- y[-zero.idx]
        w <- w[-zero.idx, -zero.idx]
    }    
    y <- y - mean(y)
    ylag <- as.numeric(w %*% y)
    sub <- paste0("MC = ", round(mc(y, w),3))
    ggplot(data.frame(y = y,
                      ylag = ylag)) +
    geom_hline(yintercept = mean(ylag),
               lty = 3) +
    geom_vline(xintercept = mean(y),
               lty = 3) +
    geom_point(
        pch = 20,
        colour = col,
        size = size,
        alpha = alpha,
        aes(x = y,
            y = ylag)
    ) +
        geom_smooth(aes(x = y,
                        y = ylag),
                method = "lm",
                lwd = lwd,
                col = "black",
                se = FALSE) +
    labs(x = xlab,
         y = ylab,
         subtitle = sub) +
    theme_classic() 
}

#' Local Moran's I
#'
#' @export
#' @description A local indicator of spatial association (lisa).
#'
#' @param x Numeric vector
#' @param w An n x n spatial connectivity matrix. See \link[geostan]{shape2mat}. This will automatically be row-standardized!
#' @param type Return the type of association also (High-High, Low-Low, High-Low, and Low-High)? Defaults to \code{FALSE}.
#'
#' @details The values will be standardized with \code{scale(x)} first and \code{w} will be row-standardized. Then the LISA is the product of a value with its mean surrounding value. These are for exploratory analysis and model diagnostics only, not ``cluster detection.'' Values greater than 2 or less than -2 are generally of interest but there is not a sound basis for applying normal distribution theory here.
#'
#' An above-average value (i.e. positive z-value) with positive mean spatial lag is of type "High-High"; a low value surrounded by high values is of type "Low-High", and so on.
#' 
#' @return If \code{type = FALSE} a numeric vector of lisa values for exploratory analysis of local spatial autocorrelation. If \code{type = TRUE}, a \code{data.frame} with columns \code{zi} (the lisa value) and \code{type}.
#'
#' @seealso \link[geostan]{moran_plot} \link[geostan]{mc} \link[geostan]{aple}
#' 
lisa <- function(x, w, type = FALSE) {
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
        message(length(zero.idx), " observations with no neighbors found. They will be converted to NA.")
    }
    w <- w / rowSums(w)
    z <- scale(x)
    lag <- as.numeric(w %*% z)
    zi <- as.numeric(z * lag)
    zi[zero.idx] <- NA
    if (!type) return (zi)
    type <- ifelse(z > 0 & lag > 0, "HH",
         ifelse(z < 0 & lag > 0, "LH",
         ifelse(z < 0 & lag < 0, "LL",
         ifelse(z>0 & lag < 0, "HL", NA
                ))))
    return (data.frame(zi = zi, type = type))
}

#' Spatial data diagnostics
#'
#' @description Visual diagnostics for areal data and model residuals
#' 
#' @param y Either a numeric vector or a \code{geostan_fit} model object (as returned from a call to one of the \code{geostan::stan_*} functions).
#' @param shape An object of class \code{sf} or another spatial object coercible to \code{sf} with \code{sf::st_as_sf} such as \code{SpatialPolygonsDataFrame}.
#' @param name The name to use on the plot labels; default to "y" or, if \code{y} is a \code{geostan_fit} object, to "Residuals".
#' @param w An optional spatial connectivity matrix; if not provided, then a row-standardized adjacency matrix will be created using \code{shape2mat(shape, "W")}.
#' @param plot If \code{FALSE}, return a list of \code{gg} plots.
#'
#' @return A grid of spatial diagnostic plots including a Moran plot plus a map and histogram of \code{y}. If a fitted \code{geostan} model is provided, model residuals are plotted and mapped (i.e. \code{y = resid(fit)$mean}).
#'
#' @seealso \link[geostan]{me_diag} \link[geostan]{mc} \link[geostan]{moran_plot} \link[geostan]{aple}
#' 
#' @export
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(sf)
#' data(ohio)
#' sp_diag(ohio$gop_growth, ohio)
#' fit <- stan_glm(gop_growth ~ 1, data = ohio, chains = 3, iter = 1e3)
#' sp_diag(fit, ohio)
#' }
#' 
sp_diag <- function(y,
                   shape,
                   name = "y",
                   w = shape2mat(shape, "W"),
                   plot = TRUE
                   ) {
    if (inherits(y, "geostan_fit")) {
        y <- resid(y)$mean
        name <- "Residuals"
        }
    if (!inherits(shape, "sf")) shape <- sf::st_as_sf(shape)
    hist <- ggplot() +
        geom_histogram(aes(y),
                       fill = "gray20",
                       col = "gray90")  +
        theme_classic() +
        labs(x = name)
    map.y <- ggplot(shape) +
        geom_sf(aes(fill = y),
                lwd = 0.05,
                col = "gray20") +
        scale_fill_gradient2(name = name) +
        theme_void()
    g.mc <- moran_plot(y, w, xlab = name)
    if (plot) {
        gridExtra::grid.arrange(hist, g.mc, map.y, ncol = 3)
    } else return (list(hist, map.y, global))
 }


#' Data model diagnostics
#'
#' @description Visual diagnostics for spatial measurement error models. 
#' 
#' @param fit A \code{geostan_fit} model object as returned from a call to one of the \code{geostan::stan_*} functions.
#' @param varname Name of the modeled variable (a character string, as it appears in the model formula).
#' @param shape An object of class \code{sf} or another spatial object coercible to \code{sf} with \code{sf::st_as_sf} such as \code{SpatialPolygonsDataFrame}.
#' @param w An optional spatial connectivity matrix; if not provided, then a row-standardized adjacency matrix will be created using \code{shape2mat(shape, "W")}.
#' @param probs Lower and upper quantiles of the credible interval to plot. 
#' @param plot If \code{FALSE}, return a \code{data.frame} with the raw data values and posterior summary of the modeled variable.
#' @param size Size of points and lines, passed to \code{geom_pointrange}.
#' 
#' @return A grid of spatial diagnostic plots for measurement error (i.e. data) models comparing the raw observations to the posterior distribution of the true values (given the specified observations, standard errors, and model). The Moran scatter plot and map depict the difference between the posterior means and the raw observations (i.e. shrinkage). 
#'
#' @seealso \link[geostan]{sp_diag} \link[geostan]{moran_plot} \link[geostan]{mc} \link[geostan]{aple}
#' @export
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(sf)
#' data(ohio)
#' C <- shape2mat(ohio)
#' ME <- list(se = data.frame(unemployment.acs = ohio$unemployment.acs.se),
#'            spatial = TRUE
#'          )
#' fit <- stan_glm(gop_growth ~ unemployment.acs,
#'                 ME = ME,
#'                 C = C,
#'                 data = ohio,
#'                 chains = 3,
#'                 iter = 1e3,
#'                 prior_only = TRUE,
#'                 refresh = 0
#'                 )
#' me_diag(fit, "unemployment.acs", ohio)
#' }
#' 
me_diag <- function(fit,
                    varname,                    
                    shape,
                    w,
                    probs = c(0.025, 0.975),
                    plot = TRUE,
                    size = 0.25
                    ) {
    if (!varname %in% colnames(fit$data)) stop("varname is not found in colnames(fit$data). Provide the name of the variable as it appears in the model formula")
    if (length(varname) != 1) stop("Provide the name of the variable as it appears in the model formula.")
    if (!inherits(shape, "sf")) shape <- sf::st_as_sf(shape)
    x.raw <- as.numeric(fit$data[,varname])
    probs = sort(probs)
    width = paste0(100 * (probs[2] - probs[1]), "%")     
    p.lwr <- paste0(100 * probs[1], "%")
    p.upr <- paste0(100 * probs[2], "%")
    x <- as.matrix(fit)
    x.samples <- x[,grep(paste0("^x_", varname), colnames(x))]
    x.mu <- apply(x.samples, 2, mean)
    x.ci <- apply(x.samples, 2, quantile, probs = probs)
    x.lwr <- x.ci[p.lwr, ]
    x.upr <- x.ci[p.upr, ] 
    df <- data.frame(
        x.raw = x.raw,
        x.mu = x.mu,
        x.lwr = x.lwr,
        x.upr = x.upr
    )
    if (plot) {
        xlbl <- paste(varname, "(raw data)")
        g.points <-  ggplot(df) +
            geom_abline(
                slope = 1,
                intercept = 0,
                lty = 2
            ) +
            geom_pointrange(
                aes(x = x.raw,
                    y = x.mu,
                    ymin = x.lwr,
                    ymax = x.upr
                    ),
                size = size
            ) +
            labs(
                x = xlbl,
                y = paste("Posterior Mean and", width, "C.I.")
            ) +                
            theme_classic()
        if (!missing(w)) {
            g.mc <- moran_plot(x.raw - x.mu, w) +
                labs(x = paste("Posterior mean minus raw data"))
        } else {
            w <- shape2mat(shape, style = "W")
            g.mc <- moran_plot(x.raw - x.mu, w) +
                labs(x = paste("Posterior mean minus raw data"))            
        }
        map.delta <- ggplot(shape) +
            geom_sf(aes(fill = df$x.mu - df$x.raw),
                    lwd = 0.05,
                    col = "gray20") +
            scale_fill_gradient2(name = paste("Posterior mean minus raw data")) +
            theme_void() +
            theme(
                legend.position = "bottom"
                ) 
        return (gridExtra::grid.arrange(g.points, g.mc, map.delta, ncol = 3))
    } else return(df)
}

#' Extract eigenfunctions of a connectivity matrix for spatial filtering
#'
#' @export
#' @param C A binary spatial weights matrix. See \link[geostan]{shape2mat} or \link[spdep]{nb2mat}.
#' @param nsa Logical. Default of \code{nsa = FALSE} excludes eigenvectors capturing negative spatial autocorrelation.
#'  Setting \code{nsa = TRUE} will result in a candidate set of EVs that contains eigenvectors representing positive and negative SA.
#' @param threshold Defaults to \code{threshold=0.2} to exclude eigenvectors representing spatial autocorrelation levels that are less than \code{threshold} times the maximum possible Moran coefficient achievable for the given spatial connectivity matrix. If \code{theshold = 0}, the eigenvector of constants (with eigenvalue of zero) will be dropped automatically.
#' @param values Should eigenvalues be returned also? Defaults to \code{FALSE}.
#' @details Returns a set of EVs limited to those with |MC| > \code{threshold} if \code{nsa = TRUE} or MC > \code{threshold} if \code{nsa = FALSE}, along with corresponding eigenvalues (optionally). Given a spatial connectivity matrix C, the function returns eigenvectors from a transformed spatial weights matrix.
#' @source
#'
#' Daniel Griffith and Yongwan Chun. 2014. "Spatial Autocorrelation and Spatial Filtering." in M. M. Fischer and P. Nijkamp (eds.), \emph{Handbook of Regional Science.} Springer.
#'
#' @return A \code{data.frame} of eigenvectors for spatial filtering. If \code{values=TRUE} then a named list is returned with elements \code{eigenvectors} and \code{eigenvalues}.
#'
#' @seealso \link[geostan]{stan_esf} \link[geostan]{mc}
#' @examples
#' 
#' library(ggplot2)
#' library(sf)
#' data(ohio)
#' C <- shape2mat(ohio, style = "B")
#' EV <- make_EV(C)
#' ggplot(ohio) +
#'   geom_sf(aes(fill = EV[,1])) +
#'   scale_fill_gradient2()
#' 
make_EV <- function(C, nsa = FALSE, threshold = 0.2, values = FALSE) {
  if (!isSymmetric(C)) C <- (t(C) + C) / 2
  N <- nrow(C)
  M <- diag(N) - matrix(1, N, N)/N
  MCM <- M %*% C %*% M
  eigens <- eigen(MCM, symmetric = TRUE)
  if(nsa) {
    idx = abs(eigens$values/eigens$values[1]) >= threshold
  } else idx <- eigens$values/eigens$values[1] >= threshold
  v <- round(eigens$values / eigens$values[1], 12)
  if (any(v == 0)) {
    rmdx <- which(v == 0)
    idx[rmdx] <- FALSE
  }
  EV <- as.data.frame(eigens$vectors[ , idx] )
  colnames(EV) <- paste0("EV", 1:ncol(EV))
  if (values) {
    lambda <- eigens$values[idx]
    return(list(eigenvectors = EV, eigenvalues = lambda))
  } else
    return(EV)
}

#' Create spatial and space-time connectivity matrices
#'
#' @export
#' @import spdep
#' @description A wrapper function for a string of \link{spdep} (and other) functions required to convert spatail objects to spatial or spatio-temporal connectivity matrices.
#' @param shape An object of class \code{sf}, \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}.
#' @param style What kind of coding scheme should be used to create the spatial connectivity matrix? Defaults to "B" for binary; use "W" for row-standardized weights; "C" for globally standardized and "S" for the Tiefelsdorf et al.'s (1999) variance-stabilizing scheme. This is passed internally to \link[spdep]{nb2mat}.
#' @param t Number of time periods. Currently only the binary coding scheme is available for space-time connectivity matrices.
#' @param st.type For space-time data, what type of space-time connectivity structure should be used? Options are "lag" for the lagged specification and "contemp" (the default) for contemporaneous specification.
#' @param zero.policy Are regions with zero neighbors allowed? Default \code{zero.policy = TRUE} (allowing regions to have zero neighbors). Also passed to \link[spdep]{nb2mat}.
#' @param queen Passed to \link[spdep]{poly2nb} to set the contiguity condition. Defaults to \code{TRUE} so that a single shared boundary point between polygons is sufficient for them to be considered neighbors.
#' @param snap Passed to \link[spdep]{poly2nb}; "boundary points less than ‘snap’ distance apart are considered to indicate contiguity."
#' 
#' @return A spatial connectivity matrix
#'
#' @seealso \link[geostan]{edges}
#' 
#' @details
#'
#' Haining and Li (Ch. 4) provide a helpful discussion of spatial connectivity matricies (Ch. 4) and spatio-temporal connectivity matricies (Ch. 15).
#'
#' The `lagged' space-time structure connects each observation to its own past (one period lagged) value and the past value of its neighbors. The `contemporaneous' specification links each observation to its neighbors and to its own in situ past (one period lagged) value (Griffith 2012, p. 23).
#' @source
#'
#' Griffith, D. A. (2012). Space, time, and space-time eigenvector filter specifications that account for autocorrelation. Estadística Espanola, 54(177), 7-34.
#' 
#' Griffith, D. A., Chun, Y., Li, B. (2020). Spatial Regression Analysis Using Eigenvector Spatial Filtering. Academic Press, Ch. 8.
#'
#' Haining, R. P., & Li, G. (2020). Regression Modelling Wih Spatial and Spatial-Temporal Data: A Bayesian Approach. CRC Press.
#' 
#'
#' @examples
#' data(ohio)
#' C <- shape2mat(ohio, "B")
#' W <- shape2mat(ohio, "W")
#'
#' ## for space-time data
#' ## if you have multiple years with same neighbors
#' ## provide the geography (for a single year!) and number of years \code{t}
#' ## defaults to the contemporaneous connectivity structure
#' Cst <- shape2mat(ohio, t = 5)
#' 
shape2mat <- function(shape, style = c("B", "W", "C", "S"), t = 1, st.type = "contemp", zero.policy = TRUE, queen = TRUE, snap = sqrt(.Machine$double.eps)) {
  style <- match.arg(style)
  shape_class <- class(shape)
  if (!any(c("sf", "SpatialPolygonsDataFrame", "SpatialPolygons") %in% shape_class)) stop("Shape must be of class SpatialPolygons, SpatialPolygonsDataFrame, or sf (simple features).")
      w <- spdep::nb2mat(spdep::poly2nb(shape, queen = queen, snap = snap), style = style, zero.policy = zero.policy)
  attributes(w)$dimnames <- NULL
  if (t > 1) { 
      if (style != "B") stop ("Only the binary coding scheme (style = 'B') has been implemented for space-time matrices.")
      ## binary temporal connectivity matrix
      s <- nrow(w)
      Ct <- matrix(0, nrow = t, ncol = t)
      for (i in 2:t) Ct[i, i-1] <- Ct[i-1, i] <- 1
      if (st.type == "lag") w = kronecker(Ct, (w + diag(s)))
      if (st.type == "contemp") {
          ## create identify matrices for space and time
          It <- diag(1, nrow = t)
          Is <- diag(1, nrow = s)
          w <- kronecker(It, w) + kronecker(Ct, Is)
      }
      }
  return(w)
}

#' Student t family
#'
#' @export
#' @description create a family object for the Student t likelihood
#' @return An object of class \code{family}
#'
student_t <- function() {
  family <- list(family = "student_t", link = 'identity')
  class(family) <- "family"
  return(family)
}

#' WAIC
#'
#' @description Widely Application Information Criteria (WAIC) for model evalution
#' @export
#' @param fit An \code{geostan_fit} object or any Stan model with a parameter named "log_lik", the pointwise log predictive likelihood
#' @param pointwise Logical, should a vector of values for each observation be returned? Default is \code{FALSE}.
#' @param digits Defaults to 2. Round results to this many digits.
#' @return A vector of length 3 with \code{WAIC}, a rough measure of the effective number of parameters estimated by the model \code{Eff_pars}, and log predictive density (\code{Lpd}). If \code{pointwise = TRUE}, results are returned in a \code{data.frame}.
#' @seealso \link{loo}
#'
#' @examples
#' 
#' \dontrun{
#' data(ohio)
#' fit <- stan_glm(gop_growth ~ 1, data = ohio, chains = 3, iter = 1500)
#' waic(fit)
#' }
#' 
waic <- function(fit, pointwise = FALSE, digits = 2) {
  ll <- as.matrix(fit, pars = "log_lik")
  nsamples <- nrow(ll)
  lpd <- apply(ll, 2, log_sum_exp) - log(nsamples)
  p_waic <- apply(ll, 2, var)
  waic <- -2 * (lpd - p_waic)
  if(pointwise) return(data.frame(waic = waic, eff_pars = p_waic, lpd = lpd))
  res <- c(WAIC = sum(waic), Eff_pars = sum(p_waic), Lpd = sum(lpd))
  return(round(res, digits))
}

#' Expected value of the residual Moran coefficient.
#'
#' @description Expected value for the Moran coefficient of model residuals under the null hypothesis of no spatial autocorrelation.
#' @export
#' @param X model matrix, including column of ones.
#' @param C Connectivity matrix.
#' @source
#'  Chun, Yongwan and Griffith, Daniel A. (2013). Spatial statistics and geostatistics. Sage, p. 18.
#' @return Returns a numeric value.
#'
expected_mc <- function(X, C) {
    n = nrow(X)
    k = ncol(X)
    under <- (n-k) * sum(rowSums(C))
    mc = -n * sum(diag( solve(t(X) %*% X) %*% t(X) %*% C %*% X )) / under
    return(as.numeric(mc))
}

#' Expected dimensions of an eigenvector spatial filter
#'
#' @description Provides an informed guess for the number of eigenvectors required to remove spatial autocorrelation from a regression.
#' For \link[geostan]{stan_esf} the result can be
#' used to set the hyper parameter \code{p0}, controlling the hyper-prior scale parameter for the global shrinkage parameter in the regularized horseshoe prior.
#' @export
#' @importFrom stats model.matrix residuals lm
#' @param formula Model formula.
#' @param data The data used to fit the model; must be coercible to a dataframe for use in \code{model.matrix}.
#' @param C An N x N binary connectivity matrix.
#' @return Returns a numeric value representing the expected number of eigenvectors required to estimate a spatial filter (i.e. number of non-zero or 'large' coefficients).
#' 
#' @details Following Chun et al. (2016), the expected number of eigenvectors required to remove residual spatial autocorrelation from a model
#'  is an increasing function of the degree of spatial autocorrelation in the outcome variable and the number of links in the connectivity matrix.
#'
#' @seealso \link[geostan]{stan_esf}
#' 
#' @source
#'
#' Chun, Yongwan, Griffith, Daniel A., Lee, Mongyeon, and Sinha, Parmanand (2016). "Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters." Journal of Geographical Systems 18(1): 67-85.
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. Spatial Statistics. \link{https://doi.org/10.1016/j.spasta.2020.100450}
#' 
exp_pars <- function(formula, data, C) {
  nlinks <- length(which(C != 0))
  N <- nrow(C)
    if (any(!C %in% c(0, 1))) {
        C <- apply(C, 2, function(i) ifelse(i != 0, 1, 0))
      }
  M <- diag(N) - matrix(1, N, N)/N
  MCM <- M %*% C %*% M
  eigens <- eigen(MCM, symmetric = TRUE)
  npos <- sum(eigens$values > 0)
  sa <- mc(residuals(lm(formula, data = data)), C)
  X <- model.matrix(formula, data)
  E_sa <- expected_mc(X, C)
  Sigma_sa <- sqrt( 2 / nlinks )
  z_sa <- (sa - E_sa) / Sigma_sa
  if (z_sa < -.59) {
    z_sa = -.59
  warning("The moran coefficient indicates very strong negative spatial autocorrelation, which this formula for obtaining the expected no. of eigenvectors was not designed for.")
  }
  a <- (6.1808 * (z_sa + .6)^.1742) / npos^.1298
  b <- 3.3534 / (z_sa + .6)^.1742
  denom <- 1 + exp(2.148 - a + b)
  candidates <- round(npos / denom)
  return(candidates)
}

#' Edge list
#'
#' @description Creates a list of connected nodes following the graph representation of a spatial connectivity matrix.
#' @export
#' @param w A connectivity matrix where connection between two nodes is indicated by non-zero entries.
#' @return Returns a \code{data.frame} with two columns representing connected pairs of nodes; only unique pairs of nodes are included.
#'
#' @details This is used internally for  \link[geostan]{stan_icar} and it is also helpful for creating the scaling factor for BYM2 models fit with \code{stan_icar}.
#'
#' @seealso \link[geostan]{shape2mat}
#' @examples
#' 
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' nbs <- edges(C)
#' head(nbs)
#' 
edges <- function(w) {
  lw <- apply(w, 1, function(r) {
    which(r != 0)
  })
  all.edges <- lapply(1:length(lw), function(i) {
    nbs <- lw[[i]]
    if(length(nbs)) data.frame(node1 = i, node2 = nbs)
  })
  all.edges <- do.call("rbind", all.edges)
  edges <- all.edges[which(all.edges$node1 < all.edges$node2),]
  return(edges)
}

#' Standard error of log(x)
#'
#' @description Transform the standard error of \code{x} to standard error of \code{log(x)}.
#'
#' @param x Estimated value of the variable \code{x}
#' @param se Standard error of \code{x}
#' @param method \code{"delta"} method uses a Taylor series approximation; the default method \code{"mc"} uses a monte carlo method.
#' @param nsim Number of draws to take if \code{method = "mc"}.
#' @param bounds Lower and upper bounds for the variable, used in the monte carlo method. Must be a length-two numeric vector with lower bound greater than or equal to zero (i.e. \code{c(lower, upper)} as in default \code{bounds = c(0, Inf)}.
#' @details The delta method returns \code{x^(-1) * se}. The monte carlo method is detailed in the examples section.
#'
#' @export
#' @importFrom truncnorm rtruncnorm
#' @examples
#' 
#' data(ohio)
#' x = ohio$unemployment.acs
#' se = ohio$unemployment.acs.se
#'
#' lse1 = se_log(x, se)
#' lse2 = se_log(x, se, method = "delta")
#' plot(lse1, lse2); abline(0, 1)
#'
#' # the monte carlo method
#' x = 10
#' se = 2
#' z = rnorm(n = 30e3, mean = x,  sd = se)
#' l.z = log(z)
#' sd(l.z)
#' se_log(x, se, method = "mc")
#' se_log(x, se, method = "delta")
#' 
se_log <- function(x, se, method = c("mc", "delta"), nsim = 30e3, bounds = c(0, Inf)) {
    stopifnot(length(x) == length(se))
    method <- match.arg(method)
    if (method == "mc") {
        stopifnot(bounds[1] >= 0 &
                  bounds[1] < bounds[2])
        se.log <- NULL
        for (i in seq_along(x)) {            
            z <- truncnorm::rtruncnorm(n = nsim, mean = x[i], sd = se[i],
                                       a = bounds[1], b = bounds[2])            
            l.z <- log(z)
            se.log <- c(se.log, sd(l.z))            
        }
        return (se.log)
    }
    if (method == "delta") return (x^(-1) * se)
}

#' Prepare data for ICAR models
#'
#' @description Given a connectivity matrix, prepare data for intrinsic conditional autoregressive models in Stan.
#' 
#' @param C connectivity matrix
#' @param scale_factor n-length vector with the scale factor for each observation's respective group. If not provided by the user it will be fixed to \code{rep(1, n)}
#' 
#' @importFrom spdep poly2nb n.comp.nb
#' 
#' @return list of data to add to Stan data list:
#' \describe{
#' \item{k}{number of groups}
#' \item{group_size}{number of nodes per group}
#' \item{n_edges}{number of connections between nodes (unique pairs only)}
#' \item{node1}{first node}
#' \item{node2}{second node. (node1[i] and node2[i] form a connected pair)}
#' \item{group_idx}{indices for each observation belonging each group, ordered by group.}
#' }
#'
#' @details This is used internally to prepare data for \link[geostan]{stan_iar} models. It can also be helpful for fitting custom ICAR models outside of \code{geostan}.
#' 
#' @seealso \link[geostan]{stan_icar} \link[geostan]{edges} \link[geostan]{shape2mat}
#'
#' @examples
#' 
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' icar.data.list <- prep_icar_data(C)
#' 
#' @export
#' @importFrom spdep n.comp.nb graph2nb
prep_icar_data <- function(C, scale_factor = NULL) {
    n <- nrow(C)
    if (inherits(scale_factor, "NULL")) scale_factor <- rep(1, n)
    E <- edges(C)
    G <- list(np = nrow(C), # confrom to spdep graph structure
              from = E$node1,
              to = E$node2,
              nedges = nrow(E)
              )
    class(G) <- "Graph"
    nb2 <- spdep::n.comp.nb(spdep::graph2nb(G))
    k = nb2$nc
    group_idx = NULL
    for (j in 1:k) group_idx <- c(group_idx, which(nb2$comp.id == j))
    group_size <- NULL
    for (j in 1:k) group_size <- c(group_size, sum(nb2$comp.id == j))
    l <- list(
        k = k,
        group_size = array(group_size, dim = k),
        n_edges = nrow(E),
        node1 = E$node1,
        node2 = E$node2,
        group_idx = array(group_idx, dim = n),
        scale_factor = scale_factor,
        comp.id = nb2$comp.id
    )
    return (l)
}

