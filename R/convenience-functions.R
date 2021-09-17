#' Effective sample size 
#'
#' @description An approximate calculation for the effective sample size for spatially autocorrelated data. Only valid for approximately normally distributed data.
#' 
#' @param n Number of observations.
#' @param rho Spatial autocorrelation parameter from a simultaneous autoregressive model.
#' 
#' @return Returns effective sample size n*, a numeric value.
#'
#' @details
#'
#' Implements Equation 3 from Griffith (2005). 
#'
#' @seealso \link[geostan]{sim_sar}, \link[geostan]{aple}
#'
#' @source
#'
#' Griffith, Daniel A. (2005). Effective geographic sample size in the presence of spatial autocorrelation. Annals of the Association of American Geographers. Vol. 95(4): 740-760.
#' 
#' @examples
#'
#' n_eff(100, 0)
#' n_eff(100, 0.5)
#' n_eff(100, 0.9)
#' n_eff(100, 1)
#'
#' rho <- seq(0, 1, by = 0.01)
#' plot(rho, n_eff(100, rho), type = 'l')
#' 
#' @export
#' 
n_eff <- function(n, rho) {
    a = 1 / (1 - exp(-1.92369)) 
    b = (n-1) / n
    c = (1 - exp(-2.12373 * rho + 0.20024 * sqrt(rho)))
    n_eff <- n * (1 - a * b * c)
    return (n_eff)
}

#' Spatial autocorrelation estimator
#'
#' @export
#'
#' @md
#' 
#' @description The approximate-profile likelihood estimator for the spatial autocorrelation parameter from a simultaneous autoregressive (SAR) model. Note, the `APLE` approximation is quite unreliable when the number of observations is large.
#' 
#' @param x Numeric vector of values, length `n`. This will be standardized internally with \code{scale(x)}.
#' @param w An `n x n` row-standardized spatial connectivity matrix. See \link[geostan]{shape2mat}.
#' @param digits Number of digits to round results to.
#' 
#' @return the APLE estimate.
#'
#' @seealso \link[geostan]{mc}, \link[geostan]{moran_plot}, \link[geostan]{lisa}, \link[geostan]{sim_sar}
#'
#' @details To check reliability, the \code{APLE} can be compared to an estimate of the spatial autocorrelation parameter from an intercept-only SAR model.
#'
#' @source
#'
#' Li, Honfei and Calder, Catherine A. and Cressie, Noel (2007). Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis: 39(4): 357-375.
#' 
#' @examples
#' 
#' library(sf)
#' data(georgia)
#' w <- shape2mat(georgia, "W")
#' x <- georgia$ICE
#' aple(x, w)
#'
#' @importFrom Matrix t
aple <- function(x, w, digits = 3) {
    check_sa_data(x, w)
    if(any(round(Matrix::rowSums(w), 10) != 1)) w <- row_standardize(w)
    z <- as.numeric(scale(x))
    n <- length(z)
    I <- diag(n)
    lambda <- eigen(w)$values
    w2 <- (w + Matrix::t(w)) / 2
    top <- Matrix::t(z) %*% w2 %*% z
    wl <- (Matrix::t(w) %*% w + as.numeric(Matrix::t(lambda) %*% lambda) * I / n)
    bottom <- Matrix::t(z) %*% wl %*% z
    return( round(as.numeric( top / bottom ), digits = digits) )
}

#' Simulate spatially autocorrelated data
#'
#' @description Given a spatial weights matrix and degree of autocorrelation, returns autocorrelated data. 
#'
#' @export
#'
#' @md
#' 
#' @param m The number of samples required. Defaults to \code{m=1} to return an \code{n}-length vector; if \code{m>1}, an \code{m x n} matrix is returned (i.e. each row will contain a sample of correlated values).
#' @param mu An \code{n}-length vector of mean values. Defaults to a vector of zeros with length equal to \code{nrow(w)}.
#' @param w Row-standardized \code{n x n} spatial weights matrix.
#' @param rho Spatial autocorrelation parameter in the range (-1, 1). Typically a scalar value; otherwise an n-length numeric vector.
#' @param sigma Scale parameter (standard deviation). Defaults to \code{sigma = 1}. Typically a scalar value; otherwise an n-length numeric vector.
#' @param ... further arguments passed to \code{MASS::mvrnorm}.
#' 
#' @return
#'
#' If \code{m = 1} a vector of the same length as \code{mu}, otherwise an \code{m x length(mu)} matrix with one sample in each row.
#'
#' @details Calls \code{MASS::mvrnorm} internally to draw from the multivariate normal distribution. The covariance matrix is specified following the simultaneous autoregressive (SAR) model. 
#'
#' @seealso \code{\link[geostan]{aple}}, \code{\link[geostan]{mc}}, \code{\link[geostan]{moran_plot}}, \code{\link[geostan]{lisa}}, \code{\link[geostan]{shape2mat}}
#' 
#' @examples
#' 
#' data(georgia)
#' w <- shape2mat(georgia, "W")
#' x <- sim_sar(w=w, rho=.5)
#' aple(x, w)
#'
#' @importFrom MASS mvrnorm
#' 
sim_sar <- function(m = 1, mu = rep(0, nrow(w)), w, rho, sigma = 1, ...) {
    check_sa_data(mu, w)
    stopifnot(all(round(Matrix::rowSums(w), 10) == 1))
    N <- nrow(w)
    if (!inherits(rho, "numeric") || !length(rho) %in% c(1, N) || any(rho >= 1 || rho <= -1)) stop("rho must be numeric value within range (-1, 1), or a k-length numeric vector where K=nrow(W).")
    if (!inherits(sigma, "numeric") || !length(sigma) %in% c(1, N) || !all(sigma > 0)) stop("sigma must be a positive numeric value, or k-length numeric vector, with K=nrow(W).")
    I <- diag(N)
    S <- sigma^2 * solve( (I - rho * Matrix::t(w)) %*% (I - rho * w) )
##    S <- crossprod(solve(I - rho * w)) * sigma^2
    x <- MASS::mvrnorm(n = m, mu = mu, Sigma = S, ...)
    return(x)
}

#' The Moran coefficient
#'
#' @description The Moran coefficient, a measure of spatial autocorrelation (also known as Global Moran's I)
#' 
#' @export
#' 
#' @param x Numeric vector of input values, length n.
#' @param w An n x n spatial connectivity matrix. See \link[geostan]{shape2mat}. 
#' @param digits Number of digits to round results to.
#' @param warn If `FALSE`, no warning will be printed to inform you when observations with zero neighbors have been dropped. 
#' 
#' @return The Moran coefficient, a numeric value.
#'
#' @details If any observations with no neighbors are found (i.e. \code{any(Matrix::rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#'
#' @seealso \link[geostan]{moran_plot}, \link[geostan]{lisa}, \link[geostan]{aple}
#' 
#' @examples
#' 
#' library(sf)
#' data(georgia)
#' w <- shape2mat(georgia, style = "W")
#' x <- georgia$ICE
#' mc(x, w)
#'
#' @source
#'
#' Chun, Yongwan, and Daniel A. Griffith. Spatial statistics and geostatistics: theory and applications for geographic information science and technology. Sage, 2013.
#' 
#' Cliff, Andrew David, and J. Keith Ord. Spatial processes: models & applications. Taylor & Francis, 1981.
#' @importFrom Matrix rowSums
#' 
mc <- function(x, w, digits = 3, warn = TRUE) {
    check_sa_data(x, w)
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        if (warn) message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data.")
        x <- x[-zero.idx]
        w <- w[-zero.idx, -zero.idx]
    }
    xbar <- mean(x)
    z <- x - xbar
    ztilde <- as.numeric(w %*% z)
    A <- sum(Matrix::rowSums(w))
    n <- length(x)
    mc <- as.numeric( n/A * (z %*% ztilde) / (z %*% z))
    return(round(mc, digits = digits))
}

#' Moran plot
#'
#' @description Plots a set of values against their spatially lagged values and gives the Moran coefficient as a measure of spatial autocorrelation.
#' 
#' @export
#' 
#' @import ggplot2
#' @param x A numeric vector of length n.
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
#' If any observations with no neighbors are found (i.e. \code{any(Matrix::rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#' 
#' @return Returns a \code{gg} plot, a scatter plot with \code{x} on the horizontal and its spatially lagged values on the vertical axis (i.e. a Moran plot).
#'
#' @seealso \link[geostan]{mc}, \link[geostan]{lisa}, \link[geostan]{aple}
#'
#' @source
#'
#' Anselin, Luc. "Local indicators of spatial association—LISA." Geographical analysis 27, no. 2 (1995): 93-115.
#' 
#' @examples
#' 
#' data(georgia)
#' x <- georgia$ICE
#' w <- shape2mat(georgia, "W")
#' moran_plot(x, w)
#'
#' @import ggplot2
#' @importFrom signs signs
#' @importFrom Matrix rowSums
moran_plot <- function(x, w, xlab = "x (centered)", ylab = "Spatial Lag", pch = 20, col = "darkred", size = 2, alpha = 1, lwd = 0.5) {
    check_sa_data(x, w)
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data.")
        x <- x[-zero.idx]
        w <- w[-zero.idx, -zero.idx]
    }    
    x <- x - mean(x)
    xlag <- as.numeric(w %*% x)
    sub <- paste0("MC = ", round(mc(x, w),3))
    ggplot(data.frame(x = x,
                      xlag = xlag),
           aes(x = x, y = xlag)
           ) +
        geom_hline(yintercept = mean(xlag),
                   lty = 3) +
        geom_vline(xintercept = mean(x),
                   lty = 3) +
        geom_point(
            pch = 20,
            colour = col,
            size = size,
            alpha = alpha,
            aes(x = x,
                y = xlag)
        ) +
        geom_smooth(aes(x = x,
                        y = xlag),
                    method = "lm",
                    formula = "y ~ x",
                    lwd = lwd,
                    col = "black",
                    se = FALSE) +
        scale_x_continuous(labels = signs::signs) +
        scale_y_continuous(labels = signs::signs) +
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
#' @param x Numeric vector of length `n`.
#' @param w An `n x n` spatial connectivity matrix. See \link[geostan]{shape2mat}. If \code{w} is not row standardized (\code{all(Matrix::rowSums(w) == 1)}), it will automatically be row-standardized.
#' @param type Return the type of association also (High-High, Low-Low, High-Low, and Low-High)? Defaults to \code{FALSE}.
#'
#' @details
#'
#' The values will be standardized with \code{z = scale(x)} first and \code{w} will be row-standardized if needed. The LISA values are the product of each \code{z} value with their respective mean surrounding value \code{lagz = w \%*\% z}; \code{lisa = z * lagz}. These are for exploratory analysis and model diagnostics. The function uses Equation 7 from Anselin (1995).
#'
#' An above-average value (i.e. positive z-value) with positive mean spatial lag indicates local positive spatial autocorrelation and is designated type "High-High"; a low value surrounded by high values indicates negative spatial autocorrelation and is designated type "Low-High", and so on.
#' 
#' @return If \code{type = FALSE} a numeric vector of lisa values for exploratory analysis of local spatial autocorrelation. If \code{type = TRUE}, a \code{data.frame} with columns \code{Li} (the lisa value) and \code{type}.
#'
#' @seealso \code{\link[geostan]{moran_plot}}, \code{\link[geostan]{mc}}, \code{\link[geostan]{aple}}
#'
#' @source
#'
#' Anselin, Luc. "Local indicators of spatial association—LISA." Geographical Analysis 27, no. 2 (1995): 93-115.
#'
#' @examples
#' 
#' library(ggplot2)
#' library(sf)
#' 
#' data(georgia)
#' w <- shape2mat(georgia, "W")
#' x <- georgia$ICE
#' li = lisa(x, w)
#' head(li)
#'
#' ggplot(georgia, aes(fill = li$Li)) +
#'   geom_sf() +
#'   scale_fill_gradient2()
#' @importFrom Matrix rowSums
lisa <- function(x, w, type = TRUE) {
    check_sa_data(x, w)
    if (!all(round(Matrix::rowSums(w), 10) %in% c(0, 1))) w <- row_standardize(w)
    z <- scale(x)
    lag <- as.numeric(w %*% z)
    zi <- as.numeric(z * lag)
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        zi[zero.idx] <- NA        
        message(length(zero.idx), " observations with no neighbors found. Their LISA values are being set to NA.")
    }     
    if (!type) return (zi)
    type <- ifelse(z > 0 & lag > 0, "HH",
         ifelse(z < 0 & lag > 0, "LH",
         ifelse(z < 0 & lag < 0, "LL",
         ifelse(z > 0 & lag < 0, "HL", NA
                ))))
    return (data.frame(Li = zi, type = type))
}

#' Spatial data diagnostics
#'
#' @description Visual diagnostics for areal data and model residuals
#' 
#' @param y A numeric vector.
#' @param shape An object of class \code{sf} or another spatial object coercible to \code{sf} with \code{sf::st_as_sf} such as \code{SpatialPolygonsDataFrame}.
#' @param name The name to use on the plot labels; default to "y" or, if \code{y} is a \code{geostan_fit} object, to "Residuals".
#' @param plot If \code{FALSE}, return a list of \code{gg} plots.
#' @param style Style of connectivity matrix; if `w` is not provided, `style` is passed to \code{\link[geostan]{shape2mat}} and defaults to "W" for row-standardized.
#' @param w An optional spatial connectivity matrix; if not provided, one will be created using \code{\link[geostan]{shape2mat}}.
#' 
#' @param binwidth A function with a single argument that will be passed to the `binwidth` argument in \code{\link[ggplot2]{geom_histogram}}. The default is to set the width of bins to `0.5 * sd(x)`.
#'
#' @param ... Additional arguments passed to \code{\link[geostan]{residuals.geostan_fit}}. For binomial and Poisson models, this includes the option to view the outcome variable as a rate (the default) rather than a count; for \code{\link[geostan]{stan_car}} models with auto-Gaussian likelihood (`fit$family$family = "auto_gaussian"), the residuals will be detrended by default, but this can be changed using `detrend = FALSE`.
#' 
#' @return A grid of spatial diagnostic plots. When provided with a numeric vector, this function plots a histogram, Moran scatter plot, and map. When provided with a fitted `geostan` model, the function returns a histogram of residuals, a histogram of Moran coefficient values calculated from the joint posterior distribution of the residuals, and a map of the mean posterior residuals (means of the marginal distributions).
#'
#' If `plot = TRUE`, the `ggplots` are drawn using \code{\link[gridExtra]{grid.arrange}}; otherwise, they are returned in a list. For the `geostan_fit` method, the underlying data for the Moran coefficient will also be returned if `plot = FALSE`.
#'
#' @seealso \code{\link[geostan]{me_diag}}, \code{\link[geostan]{mc}}, \code{\link[geostan]{moran_plot}}, \code{\link[geostan]{aple}}
#' 
#'
#' @examples
#' data(georgia)
#' sp_diag(georgia$college, georgia)
#'
#' bin_fn <- function(y) mad(y)
#' sp_diag(georgia$college, georgia, binwidth = bin_fn)
#' 
#' \dontrun{
#' fit <- stan_glm(log(rate.male) ~ 1, data = georgia)
#' sp_diag(fit, georgia)
#'
#' cp <- prep_car_data(shape2mat(georgia))
#' fit2 <- stan_car(log(rate.male) ~ 1, data = georgia, car_parts = cp)
#' sp_diag(fit2, georgia)
#' sp_diag(fit2, georgia, detrend = FALSE)
#' }
#'
#' @importFrom stats sd
#' @export
#' @md
sp_diag <- function(y,
                   shape,
                   name = "y",
                   plot = TRUE,
                   style = c("W", "B"),
                   w = shape2mat(shape, match.arg(style)),
                   binwidth = function(x) 0.5 * sd(x),
                   ...
                   ) {
    if (inherits(y, "integer")) y <- as.numeric(y)
    UseMethod("sp_diag", y)
}

#' @export
#' @md
#' @param y A fitted `geostan` model (class `geostan_fit`).
#' @method sp_diag geostan_fit
#' @rdname sp_diag
#' @importFrom sf st_as_sf
#' @importFrom gridExtra grid.arrange
#' @importFrom signs signs
#' @importFrom stats sd
#' @import ggplot2
sp_diag.geostan_fit <- function(y,
                            shape,
                            name = "Residual",
                            plot = TRUE,
                            style = c("W", "B"),
                            w = shape2mat(shape, match.arg(style)),
                            binwidth = function(x) 0.5 * stats::sd(x),
                            ...) {
    marginal_residual <- residuals(y, summary = TRUE, ...)$mean
    if (!inherits(shape, "sf")) shape <- sf::st_as_sf(shape)
    hist.y <- ggplot() +
        geom_histogram(aes(marginal_residual),                       
                       binwidth = binwidth(marginal_residual),
                       fill = "gray20",
                       col = "gray90")  +
        scale_x_continuous(labels = signs::signs) +
        theme_classic() +
        labs(x = name)
    map.y <- ggplot(shape) +
        geom_sf(aes(fill = marginal_residual),
                lwd = 0.05,
                col = "gray20") +
        scale_fill_gradient2(name = name,
                             label = signs::signs) +
        theme_void()   
    R <- residuals(y, summary = FALSE)
    R.mc <- apply(R, 1, mc, w = w)
    if (length(unique(R.mc)) == 1) {
        g.mc <- moran_plot(R[1,], w, xlab = name)
    } else {
        R.mc.mu <- mean(R.mc)
        g.mc <- ggplot() +
            geom_histogram(aes(R.mc),
                                        #     binwidth = binwidth(as.numeric(R.mc)),
                           fill = "gray20",
                           col = "gray90") +
            scale_x_continuous(labels = signs::signs) +
            theme_classic() +
            labs(
                x = "Residual MC",
                subtitle = paste0("MC (mean) = ", round(R.mc.mu, 2)))
    }
    if (plot) {
        return( gridExtra::grid.arrange(hist.y, g.mc, map.y, ncol = 3) )
    } else {        
        return (list(residual_histogram = hist.y, mc_plot = g.mc, residual_map = map.y, mc_data = R.mc))
        }
 }


#' @export
#' @method sp_diag numeric
#' @rdname sp_diag
#' @importFrom sf st_as_sf
#' @importFrom gridExtra grid.arrange
#' @importFrom signs signs
#' @importFrom stats sd
#' @import ggplot2
sp_diag.numeric <- function(y,
                            shape,
                            name = "Residual",
                            plot = TRUE,
                            style = c("W", "B"),
                            w = shape2mat(shape, match.arg(style)),
                            binwidth = function(x) 0.5 * stats::sd(x),
                            ...) {
    if (!inherits(shape, "sf")) shape <- sf::st_as_sf(shape)
    stopifnot(length(y) == nrow(shape))
    shape$y <- y
    hist.y <- ggplot() +
        geom_histogram(aes(y),
                       binwidth = binwidth(y),
                       fill = "gray20",
                       col = "gray90")  +
        scale_x_continuous(labels = signs::signs) +
        theme_classic() +
        labs(x = name)
    map.y <- ggplot(shape) +
        geom_sf(aes(fill = y),
                lwd = 0.05,
                col = "gray20") +
        scale_fill_gradient2(name = name,
                             label = signs::signs) +
        theme_void()
    g.mc <- moran_plot(y, w, xlab = name)
    if (plot) {
        return( gridExtra::grid.arrange(hist.y, g.mc, map.y, ncol = 3) )
    } else return (list(residual_histogram = hist.y, mc_plot = g.mc, residual_map = map.y))
 }


#' Data model diagnostics
#'
#' @description Visual diagnostics for spatial measurement error models. 
#' 
#' @param fit A \code{geostan_fit} model object as returned from a call to one of the \code{geostan::stan_*} functions.
#' @param varname Name of the modeled variable (a character string, as it appears in the model formula).
#' @param shape An object of class \code{sf} or another spatial object coercible to \code{sf} with \code{sf::st_as_sf}.
#' @param probs Lower and upper quantiles of the credible interval to plot. 
#' @param plot If \code{FALSE}, return a list of \code{ggplot}s and a \code{data.frame} with the raw data values alongside a posterior summary of the modeled variable.
#' @param size Size of points and lines, passed to \code{geom_pointrange}.
#' @param index Integer value; use this if you wish to identify observations with the largest `n=index` absolute Delta values; data on the top `n=index` observations ordered by absolute Delta value will be printed to the console and the plots will be labeled with the indices of the identified observations.
#' 
#' @param style Style of connectivity matrix; if `w` is not provided, `style` is passed to \code{\link[geostan]{shape2mat}} and defaults to "W" for row-standardized.
#' @param w An optional spatial connectivity matrix; if not provided, one will be created using \code{\link[geostan]{shape2mat}}.
#'
#' @param binwidth A function with a single argument that will be passed to the `binwidth` argument in \code{\link[ggplot2]{geom_histogram}}. The default is to set the width of bins to `0.5 * sd(x)`.
#' 
#' @return A grid of spatial diagnostic plots for measurement error models comparing the raw observations to the posterior distribution of the true values. Includes a point-interval plot, a histogram of Moran coefficient values for the posterior distribution of Delta values (`Delta = z - x`, where `z` are the survey estimates and `x` are the modeled true values), and a map of the posterior mean of the Delta values. 
#'
#' @seealso \code{\link[geostan]{sp_diag}}, \code{\link[geostan]{moran_plot}}, \code{\link[geostan]{mc}}, \code{\link[geostan]{aple}}
#' 
#' @source
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). ``Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure.'' *Int. J. Env. Res. and Public Health* 18 (13): 6856. DOI: 10.3390/ijerph18136856 Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#' 
#' @examples
#' \dontrun{
#' library(sf)
#' data(georgia)
#' ## binary adjacency matrix
#' A <- shape2mat(georgia, "B")
#' ## prepare data for the CAR model, using WCAR specification
#' cp <- prep_car_data(A, type = "WCAR")
#' ## provide list of data for the measurement error model
#' ME <- list(se = data.frame(ICE = georgia$ICE.se),
#'            spatial = TRUE,
#'            car_parts = cp
#'          )
#' ## sample from the prior probability model only, including the ME model
#' fit <- stan_glm(log(rate.male) ~ ICE,
#'                 ME = ME,
#'                 C = C,
#'                 data = georgia, 
#'                 prior_only = TRUE,
#'                 refresh = 0
#'                 )
#' ## see ME diagnostics
#' me_diag(fit, "ICE", georgia)
#' ## see index values for the largest (absolute) delta values
#'  ## (differences between raw estimate and the posterior mean)
#' me_diag(fit, "ICE", georgia, index = 3)
#' }
#'
#' @export
#' @md
#' @importFrom sf st_as_sf
#' @importFrom gridExtra grid.arrange
#' @importFrom signs signs
#' @import ggplot2
#' @importFrom stats sd
me_diag <- function(fit,
                    varname,                    
                    shape,
                    probs = c(0.025, 0.975),
                    plot = TRUE,
                    size = 0.25,
                    index = 0,
                    style = c("W", "B"),
                    w = shape2mat(shape, match.arg(style)),
                    binwidth = function(x) 0.5 * sd(x)
                    ) {
    stopifnot(length(varname) == 1)    
    if (!varname %in% colnames(fit$data)) stop("varname is not found in colnames(fit$data). Provide the name of the variable as it appears in the model formula")
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
    xlbl <- paste(varname, "(raw data)")
    g.points <-  ggplot(df,
                        aes(x = x.raw,
                            y = x.mu,
                            ymin = x.lwr,
                            ymax = x.upr
                            )
                        ) +
        geom_abline(
            slope = 1,
            intercept = 0,
            lty = 2
        ) +
        geom_pointrange(
            size = size
        ) +
        scale_x_continuous(labels = signs::signs) +
        scale_y_continuous(labels = signs::signs) +
        labs(
            x = xlbl,
            y = paste("Posterior Mean and", width, "C.I.")
        ) +                
        theme_classic()
    delta.mat <- t(apply(x.samples, 1, .resid, y = x.raw))
    df$Delta <- apply(delta.mat, 2, mean)
    D.mc <- apply(delta.mat, 1, mc, w = w)
    D.mc.mu <- mean(D.mc)
    g.mc <- ggplot() +
        geom_histogram(aes(D.mc),
                       binwidth = binwidth(as.numeric(D.mc)),
                       fill = "gray20",
                       col = "gray90") +
        scale_x_continuous(labels = signs::signs) +
        theme_classic() +
        labs(
            x = expression(paste('MC (', Delta, ')')),
            subtitle = paste0("Mean MC = ", round(D.mc.mu, 3)))    
    map.delta <- ggplot(shape) +
        geom_sf(aes(fill = df$Delta),
                lwd = 0.05,
                col = "gray20") +
        scale_fill_gradient2(name = bquote(hat(Delta)),
                             label = signs::signs) +
    theme_void() +
    theme(
        legend.position = "right"
    )
    if (index) {            
        message("Identifying the top ", index, " observations as ordered by their Delta values (Delta = posterior mean of x - raw x value):")
        ordered.ids <- order(abs(df$Delta), decreasing = TRUE)[1:index]
        print(df[ordered.ids, ])
        df$ID <- NA
        df$ID[ordered.ids] <- ordered.ids
        g.points <- g.points +
            geom_label(aes(label = df$ID),
                       na.rm = TRUE,
                       alpha = 0.9
                       )
        map.delta <- map.delta +
            geom_sf(aes(col = !is.na(df$ID),
                        fill = NULL),
                    fill = alpha("white", 0)) +
            scale_color_manual(
                    breaks =  c(FALSE, TRUE),
                    values = c(NA, "black"),
                    guide = "none"
            )
                    
    }
    if (plot) {
        return (gridExtra::grid.arrange(g.points, g.mc, map.delta, ncol = 3))
    } else {
        res.list <- list(point_interval = g.points, mc_plot = g.mc, delta_map = map.delta, delta_data = df, mc_data = D.mc)
        return(res.list)
    }
}

#' Extract eigenfunctions of a connectivity matrix for spatial filtering
#'
#' @export
#' @param C A binary spatial weights matrix. See \code{\link[geostan]{shape2mat}}.
#' @param nsa Logical. Default of \code{nsa = FALSE} excludes eigenvectors capturing negative spatial autocorrelation.
#'  Setting \code{nsa = TRUE} will result in a candidate set of EVs that contains eigenvectors representing positive and negative SA.
#' @param threshold Defaults to \code{threshold=0.2} to exclude eigenvectors representing spatial autocorrelation levels that are less than \code{threshold} times the maximum possible Moran coefficient achievable for the given spatial connectivity matrix. If \code{theshold = 0}, all eigenvectors will be returned (however, the eigenvector of constants (with eigenvalue of zero) will be dropped automatically).
#' @param values Should eigenvalues be returned also? Defaults to \code{FALSE}.
#' @details Returns a set of eigenvectors related to the Moran coefficient (MC), limited to those eigenvectors with |MC| > \code{threshold} if \code{nsa = TRUE} or MC > \code{threshold} if \code{nsa = FALSE}, optionally with corresponding eigenvalues. 
#' @source
#'
#' Daniel Griffith and Yongwan Chun. 2014. "Spatial Autocorrelation and Spatial Filtering." in M. M. Fischer and P. Nijkamp (eds.), \emph{Handbook of Regional Science.} Springer.
#'
#' @return A \code{data.frame} of eigenvectors for spatial filtering. If \code{values=TRUE} then a named list is returned with elements \code{eigenvectors} and \code{eigenvalues}.
#'
#' @seealso \link[geostan]{stan_esf}, \link[geostan]{mc}
#' @examples
#' 
#' library(ggplot2)
#' library(sf)
#' data(georgia)
#' C <- shape2mat(georgia, style = "B")
#' EV <- make_EV(C)
#' ggplot(georgia) +
#'   geom_sf(aes(fill = EV[,1])) +
#'   scale_fill_gradient2()
#'
#' \dontrun{
#' fit <- stan_esf(log(rate.male) ~ 1, data = georgia, EV = EV, C = C)
#' sp_diag(fit, georgia)
#' }
#' @importFrom Matrix isSymmetric t
make_EV <- function(C, nsa = FALSE, threshold = 0.2, values = FALSE) {
    if (!Matrix::isSymmetric(C)) {
        message("Coercing C to be symmetric with: C <- (t(C) + C) / 2")
        C <- (Matrix::t(C) + C) / 2
    }
  N <- nrow(C)
  M <- Matrix::Diagonal(N) - Matrix::Matrix(1, nrow = N, ncol = N)/N
  MCM <- M %*% C %*% M
  eigens <- eigen(MCM, symmetric = TRUE)
  if(nsa) {
    idx = abs(eigens$values/eigens$values[1]) >= threshold
  } else {
      idx <- eigens$values/eigens$values[1] >= threshold
  }
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
#' @description Creates sparse matrix representations of spatial connectivity structures
#' 
#' @param shape An object of class \code{sf}, \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}.
#' 
#' @param style What kind of coding scheme should be used to create the spatial connectivity matrix? Defaults to "B" for binary; use "W" for row-standardized weights.
#'
#' @param queen Passed to \code{\link[spdep]{poly2nb}} to set the contiguity condition. Defaults to \code{TRUE} so that a single shared boundary point (rather than a shared border/line) between polygons is sufficient for them to be considered neighbors.
#' 
#' @param snap Passed to \code{\link[spdep]{poly2nb}}; "boundary points less than ‘snap’ distance apart are considered to indicate contiguity."
#' 
#' @param t Number of time periods. Only the binary coding scheme is available for space-time connectivity matrices.
#' 
#' @param st.style For space-time data, what type of space-time connectivity structure should be used? Options are "lag" for the lagged specification and "contemp" (the default) for contemporaneous specification (see Details).
#' 
#' 
#' @return A spatial connectivity matrix
#'
#' @seealso \code{\link[geostan]{edges}}
#' 
#' @details
#'
#' Haining and Li (Ch. 4) provide a helpful discussion of spatial connectivity matrices (Ch. 4).
#'
#' The `lagged' space-time structure connects each observation to its own past (one period lagged) value and the past value of its neighbors. The `contemporaneous' specification links each observation to its neighbors and to its own in situ past (one period lagged) value (Griffith 2012, p. 23).
#' 
#' @source
#'
#' Griffith, D. A. (2012). Space, time, and space-time eigenvector filter specifications that account for autocorrelation. Estadística Espanola, 54(177), 7-34.
#'
#' Haining, R. P., & Li, G. (2020). Regression Modelling Wih Spatial and Spatial-Temporal Data: A Bayesian Approach. CRC Press.
#' 
#'
#' @examples
#' 
#' data(georgia)
#'
#' ## binary adjacency matrix
#' C <- shape2mat(georgia, "B")
#' ## row sums gives the numbers of neighbors per observation
#' Matrix::rowSums(C)
#' head(Matrix::summary(C))
#'
#' ## row-standardized matrix 
#' W <- shape2mat(georgia, "W")
#' Matrix::rowSums(W)
#' head(Matrix::summary(W))
#' 
#' ## space-time matricies 
#' ## for eigenvector space-time filtering
#' ## if you have multiple years with same neighbors,
#' ## provide the geography (for a single year!) and number of years \code{t}
#' Cst <- shape2mat(georgia, t = 5)
#' dim(Cst)
#' EVst <- make_EV(Cst)
#' dim(EVst)
#' 
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom spdep poly2nb
#' @export
#'
shape2mat <- function(shape,
                       style = c("B", "W"),
                       queen = TRUE,
                       snap = sqrt(.Machine$double.eps),
                       t = 1,
                       st.style = c("contemp", "lag"))
{
    style <- match.arg(style)
    st.style <- match.arg(st.style)
    nb <- spdep::poly2nb(shape, queen = queen, snap = snap)
    ids <- 1:nrow(shape)
    dims <- rep(length(ids), 2)
    Ni <- unlist(lapply(nb, count_neighbors))
    i <- rep(ids, times = Ni)
    j <- do.call("c", nb)
    j <- j[j>0]
    stopifnot(length(i) == length(j))
    C <- Matrix::sparseMatrix(i = i, j = j, dims = dims)
    if (style == "W") C <- row_standardize(C, warn = FALSE)
    if (t > 1) { 
        if (style != "B") stop ("Only the binary coding scheme has been implemented for space-time matrices.")
      ## binary temporal connectivity matrix
      s <- nrow(C)
      Ct <- matrix(0, nrow = t, ncol = t)
      for (i in 2:t) Ct[i, i-1] <- Ct[i-1, i] <- 1
      if (st.style == "lag") C = kronecker(Ct, (C + diag(s)))
      if (st.style == "contemp") {
          ## create identify matrices for space and time
          It <- diag(1, nrow = t)
          Is <- diag(1, nrow = s)
          C <- kronecker(It, C) + kronecker(Ct, Is)
      }
  }    
    return(C)
}


#' Used internally for shape2mat to count neighbs per element of an nb list; as in, `Ni <- unlist(lapply(nb, count_neighbors))`
#' @param z an element from spdep's nb object
#' @noRd 
count_neighbors <- function(z) {
    if (length(z) == 1 && z == 0) 0 else length(z)
}

#' Row-standardize a matrix; safe for zero row-sums.
#'
#' @param C A matrix
#' @param warn Print `msg` if `warn = TRUE`.
#' @param msg A warning message to print.
#' @return A row-standardized matrix, W (i.e., all row sums equal 1, or zero).
#' @examples
#' A <- shape2mat(georgia)
#' head(Matrix::summary(A))
#' Matrix::rowSums(A)
#' W <- row_standardize(A)
#' head(Matrix::summary(W))
#' Matrix::rowSums(W)
#' @importFrom Matrix rowSums
#' @export
row_standardize <- function(C, warn = TRUE, msg = "Row standardizing connectivity matrix") {
    stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
    if (warn) message(msg)    
    Ni <- Matrix::rowSums(C)
    Ni[Ni == 0] <- 1
    W <- C / Ni
    return (W)
}

#' Student t family
#'
#' @export
#' @description create a family object for the Student t likelihood
#' @return An object of class \code{family}
#' @examples
#' \dontrun{
#' data(georgia)
#' fit = stan_glm(log(rate.male) ~ 1, data = georgia, family = student_t())
#' }
student_t <- function() {
  family <- list(family = "student_t", link = 'identity')
  class(family) <- "family"
  return(family)
}

#' auto Gaussian family for CAR models
#'
#' @export
#' @description create a family object for the auto-Gaussian CAR specification
#' @return An object of class \code{family}
#' @seealso \code{\link[geostan]{stan_car}}
#' @examples
#' \dontrun{
#' cp = prep_car_data(shape2mat(georgia))
#' fit <- stan_car(log(rate.male) ~ 1, data = georgia, car_parts = cp, family = auto_gaussian())
#' }
auto_gaussian <- function() {
    family <- list(family = "auto_gaussian", link = "identity")
    class(family) <- "family"
    return(family)
}

#' WAIC
#'
#' @description Widely Application Information Criteria (WAIC) for model comparison
#' 
#' @param fit An \code{geostan_fit} object or any Stan model with a parameter named "log_lik", the pointwise log likelihood of the observations.
#' @param pointwise Logical (defaults to `FALSE`), should a vector of values for each observation be returned? 
#' @param digits Round results to this many digits.
#' @return A vector of length 3 with \code{WAIC}, a rough measure of the effective number of parameters estimated by the model \code{Eff_pars}, and log predictive density \code{Lpd}. If \code{pointwise = TRUE}, results are returned in a \code{data.frame}.
#' 
#' @seealso \code{\link[loo]{waic}} \code{\link[loo]{loo}}
#'
#' @examples
#' 
#' \dontrun{
#' data(georgia)
#' fit <- stan_glm(college ~ 1, data = georgia)
#' waic(fit)
#' }
#' @source
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information criterion in singular learning theory. Journal of Machine Learning Research 11, 3571-3594.
#' 
#' @export
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
#' 
#' @export
#' 
#' @param X model matrix, including column of ones.
#' @param C Connectivity matrix.
#' 
#' @source
#' 
#'  Chun, Yongwan and Griffith, Daniel A. (2013). Spatial statistics and geostatistics. Sage, p. 18.
#' 
#' @return Returns a numeric value.
expected_mc <- function(X, C) {
    C <- as.matrix(C)
    n = nrow(X)
    k = ncol(X)    
    under <- (n-k) * sum(rowSums(C))
    mc = -n * sum(diag( solve(t(X) %*% X) %*% t(X) %*% C %*% X )) / under
    return(as.numeric(mc))
}

#' Expected dimensions of an eigenvector spatial filter
#'
#' @description Provides an informed guess for the number of eigenvectors required to remove spatial autocorrelation from a regression. This is used internally for \code{\link[geostan]{stan_esf}}; the result can be used to set the prior scale parameter for the global shrinkage parameter in the regularized horseshoe prior. A smaller value of `p0` leads to a more sparse specification.
#' 
#' 
#' @param formula Model formula.
#' @param data The data used to fit the model; must be coercible to a dataframe for use in \code{model.matrix}.
#' @param C An N x N binary connectivity matrix.
#' @return Returns a numeric value representing the expected number of eigenvectors required to estimate a spatial filter (i.e. number of non-zero or 'large' coefficients).
#' 
#' @details Following Chun et al. (2016), the expected number of eigenvectors required to remove residual spatial autocorrelation from a model
#'  is an increasing function of the degree of spatial autocorrelation in the outcome variable and the number of links in the connectivity matrix.
#'
#' @seealso \code{\link[geostan]{stan_esf}}
#' 
#' @source
#'
#' Chun, Y., D. A. Griffith, M. Lee and P. Sinha (2016). Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters. *Journal of Geographical Systems*, 18(1), 67-85. \doi{10.1007/s10109-015-0225-3}.
#'
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. *Spatial Statistics*. \doi{10.1016/j.spasta.2020.100450}.
#'
#' Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In *Electronic Journal of Statistics*, 11(2):5018-5051. \doi{10.1214/17-EJS1337SI}.
#' 
#' @examples
#'
#' C <- shape2mat(georgia, "B")
#' c(p0 = exp_pars(log(rate.male) ~ college, georgia, C))
#' \dontrun{
#'  fit <- stan_esf(log(rate.male) ~ college, data = georgia, p0 = p0, iter = 1e3)
#' }
#' 
#' @importFrom stats model.matrix residuals lm
#' @export
#' @importFrom Matrix summary
exp_pars <- function(formula, data, C) {
    stopifnot(inherits(C, "matrix") | inherits(C, "Matrix"))
##    C <- as(C, "ngCMatrix")
    C <- Matrix::Matrix(C)
    nlinks <- nrow(Matrix::summary(C))
    N <- nrow(C)
    ## if (any(!C %in% c(0, 1))) {
    ##     C <- apply(C, 2, function(i) ifelse(i != 0, 1, 0))
    ## }
    M <- Matrix::Diagonal(N) - Matrix::Matrix(1, nrow = N, ncol = N)/N
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
#' 
#' @param C A connectivity matrix where connection between two nodes is indicated by non-zero entries.
#' @param unique_pairs_only By default, only unique pairs of nodes (i, j) will be included in the output.
#' 
#' @return
#' 
#' Returns a \code{data.frame} with three columns. The first two columns (\code{node1} and \code{node2}) contain the indices of connected pairs of nodes; only unique pairs of nodes are included (unless `unique_pairs_only = FALSE`). The third column (\code{weight}) contains the corresponding matrix element, \code{C[node1, node2]}.
#'
#' @details This is used internally for \code{\link[geostan]{stan_icar}} and it is also helpful for creating the scaling factor for BYM2 models fit with \code{\link[geostan]{stan_icar}}.
#'
#' @seealso \code{\link[geostan]{shape2mat}}, \code{\link[geostan]{prep_icar_data}}, \code{\link[geostan]{stan_icar}}
#' @examples
#' 
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' nbs <- edges(C)
#' head(nbs)
#'
#' ## similar to:
#' head(Matrix::summary(C))
#' head(Matrix::summary(shape2mat(georgia, "W")))
#' 
#' @importFrom Matrix summary
#' @export
edges <- function(C, unique_pairs_only = TRUE) {
    stopifnot(inherits(C, "Matrix") | inherits(C, "matrix"))
    edges <- Matrix::summary(Matrix::Matrix(C))
    names(edges)[1:2] <- c("node1", "node2")
    if ("x" %in% names(edges)){
        edges$x <- as.numeric(edges$x)
        names(edges)[3] <- "weight"
    } else {
        edges$weight <- 1
    }
    edges <- edges[order(edges$node1, edges$node2), ]
    if (unique_pairs_only) edges <- edges[which(edges$node1 < edges$node2),]
    rownames(edges) <- NULL
    class(edges) <- "data.frame"
    return(edges)
}

#' Standard error of log(x)
#'
#' @description Transform the standard error of \code{x} to standard error of \code{log(x)}.
#'
#' @param x An estimate
#' @param se Standard error of \code{x}
#' @param method The \code{"delta"} method uses a Taylor series approximation; the default method, \code{"mc"}, uses a simple monte carlo method.
#' @param nsim Number of draws to take if \code{method = "mc"}.
#' @param bounds Lower and upper bounds for the variable, used in the monte carlo method. Must be a length-two numeric vector with lower bound greater than or equal to zero (i.e. \code{c(lower, upper)} as in default \code{bounds = c(0, Inf)}.
#' @details The delta method returns \code{x^(-1) * se}. The monte carlo method is detailed in the examples section.
#'
#' @export
#' @importFrom truncnorm rtruncnorm
#' @examples
#' 
#' data(georgia)
#' x = georgia$college
#' se = georgia$college.se
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
se_log <- function(x, se, method = c("mc", "delta"), nsim = 5e3, bounds = c(0, Inf)) {
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
#' @description Given a symmetric n x n connectivity matrix, prepare data for intrinsic conditional autoregressive models in Stan. This function may be used for building custom ICAR models in Stan. This is used internally by \code{\link[geostan]{stan_icar}}.
#' 
#' @param C Connectivity matrix
#' @param scale_factor Optional vector of scale factors for each connected portion of the graph structure. If not provided by the user it will be fixed to a vector of ones. 
#' 
#' @importFrom spdep poly2nb n.comp.nb
#' 
#' @return list of data to add to Stan data list:
#' 
#' \describe{
#' \item{k}{number of groups}
#' \item{group_size}{number of nodes per group}
#' \item{n_edges}{number of connections between nodes (unique pairs only)}
#' \item{node1}{first node}
#' \item{node2}{second node. (node1[i] and node2[i] form a connected pair)}
#' \item{weight}{The element \code{C[node1, node2]}.}
#' \item{group_idx}{indices for each observation belonging each group, ordered by group.}
#' \item{m}{number of disconnected regions requiring their own intercept.}
#' \item{A}{n-by-m matrix of dummy variables for the component-specific intercepts.}
#' \item{inv_sqrt_scale_factor}{By default, this will be a k-length vector of ones. Placeholder for user-specified information. If user provided `scale_factor`, then this will be `1/sqrt(scale_factor)`.}
#' \item{comp_id}{n-length vector indicating the group membership of each observation.}
#' }
#'
#' @details
#'
#' This is used internally to prepare data for \link[geostan]{stan_icar} models. It can also be helpful for fitting custom ICAR models outside of \code{geostan}. 
#' 
#' @seealso \code{\link[geostan]{edges}}, \code{\link[geostan]{shape2mat}}, \code{\link[geostan]{stan_icar}}, \code{\link[geostan]{prep_car_data}}
#'
#' @examples
#' 
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' icar.data.list <- prep_icar_data(C)
#'
#' @source
#'
#' Besag, Julian, Jeremy York, and Annie Mollié. 1991. “Bayesian Image Restoration, with Two Applications in Spatial Statistics.” Annals of the Institute of Statistical Mathematics 43 (1): 1–20.
#' 
#' Donegan, Connor. Flexible Functions for ICAR, BYM, and BYM2 Models in Stan. Code Repository. 2021. Available online: \url{https://github.com/ConnorDonegan/Stan-IAR} (accessed Sept. 10, 2021).
#'
#' Freni-Sterrantino, Anna, Massimo Ventrucci, and Håvard Rue. 2018. “A Note on Intrinsic Conditional Autoregressive Models for Disconnected Graphs.” Spatial and Spatio-Temporal Epidemiology 26: 25–34.
#'
#' Morris, Mitzi, Katherine Wheeler-Martin, Dan Simpson, Stephen J Mooney, Andrew Gelman, and Charles DiMaggio. 2019. “Bayesian Hierarchical Spatial Models: Implementing the Besag York Mollié Model in Stan.” Spatial and Spatio-Temporal Epidemiology 31: 100301.
#'
#' Riebler, Andrea, Sigrunn H Sørbye, Daniel Simpson, and Håvard Rue. 2016. “An Intuitive Bayesian Spatial Model for Disease Mapping That Accounts for Scaling.” Statistical Methods in Medical Research 25 (4): 1145–65.
#' 
#' @export
#' @importFrom spdep n.comp.nb graph2nb
prep_icar_data <- function(C, scale_factor = NULL) {
  n <- nrow(C)
  E <- edges(C, unique_pairs_only = TRUE)
  G <- list(np = nrow(C), from = E$node1, to = E$node2, nedges = nrow(E))
  class(G) <- "Graph"
  nb2 <- spdep::n.comp.nb(spdep::graph2nb(G))
  k = nb2$nc
  if (!inherits(scale_factor, "NULL")) {
      if (length(scale_factor) != k) stop("scale_factor is of wrong length. Must have one value per fully connected graph component. See the documentation for `geostan::stan_icar` to learn how to create the scale_factor.")
      }
  if (inherits(scale_factor, "NULL")) scale_factor <- array(rep(1, k), dim = k)
  group_idx = NULL
  for (j in 1:k) group_idx <- c(group_idx, which(nb2$comp.id == j))
  group_size <- NULL
  for (j in 1:k) group_size <- c(group_size, sum(nb2$comp.id == j))
  # intercept per connected component of size > 1, if multiple.
  m <- sum(group_size > 1) - 1
  if (m) {
    GS <- group_size
    ID <- nb2$comp.id
    change.to.one <- which(GS == 1)
    ID[which(ID == change.to.one)] <- 1
    A = model.matrix(~ factor(ID))
    A <- as.matrix(A[,-1])
  } else {
    A <- model.matrix(~ 0, data.frame(a=1:n))
  }
  l <- list(k = k, 
            group_size = array(group_size, dim = k), 
            n_edges = nrow(E), 
            node1 = E$node1, 
            node2 = E$node2,
            weight = E$weight,
            group_idx = array(group_idx, dim = n), 
            m = m,
            A = A,
            inv_sqrt_scale_factor = as.array(1 / sqrt(scale_factor)),
            comp_id = nb2$comp.id)
  return(l)
}


#' Prepare data for a Stan CAR model
#'
#' 
#' @param A Binary adjacency matrix; for `style = DCAR`, provide a symmetric matrix of distances instead. The distance matrix should be sparse, meaning that most distances should be zero (usually obtained by setting some threshold distance beyond which all are zero).
#' 
#' @param style Specification for the connectivity matrix (C) and conditional variances (M); one of "WCAR", "ACAR", or "DCAR".
#'
#' @param k For `style = DCAR`, distances will be raised to the -k power (d^-k).
#'
#' @param gamma For `style = DCAR`, distances will be offset by `gamma` before raising to the `-k`th power.
#' 
#' @param lambda If TRUE, return eigenvalues required for calculating the log determinant of the precision matrix and for determining the range of permissible values of rho. These will also be printed with a message if lambda = TRUE.
#' 
#' @param cmat If `cmat = TRUE`, return the full matrix C (in sparse matrix format).
#' 
#' 
#' @details
#' The CAR model is:
#' ```
#'   Normal(Mu, Sigma), Sigma = (I - rho * C)^-1 * M * tau^2,
#' ```
#' where `I` is the identity matrix, `rho` is a spatial autocorrelation parameter, `C` is a connectivity matrix, and `M * tau^2` is a diagonal matrix with conditional variances on the diagonal. `tau^2` is a (scalar) scale parameter.
#'
#' In the WCAR specification, `C` is the row-standardized version of `A`. This means that the non-zero elements of `A` will be converted to `1/N_i` where `N_i` is the number of neighbors for the `i`th site (obtained using `Matrix::rowSums(A)`. The conditional variances (on the diagonal of `M * tau^2`), are also proportional to `1/N_i`. 
#'
#' The ACAR specification is from Cressie, Perrin and Thomas-Agnon (2005); also see Cressie and Wikle (2011, p. 188).
#'
#' The DCAR specification is inverse distance-based, and requires the user provide a (sparse) distance matrix instead of a binary adjacency matrix. (For `A`, provide a symmetric matrix of distances, not inverse distances!) Internally, non-zero elements of `A` will be converted to: `d_{ij} = (a_{ij} + gamma)^(-k)` (Cliff and Ord 1981, p. 144). Default values are `k=1` and `gamma=0`. Following Cressie (2015), these values will be standardized by the maximum `d_{ij}` value. The conditional variances will be proportional to the inverse of the row sums of the transformed distance matrix: `D_{ii} = (sum_i^N d_{ij})^(-1)`.
#'
#' For inverse-distance weighting schemes, see Cliff and Ord (1981); for distance-based CAR specifications, see Cressie (2015 \[1993\]) and Haining and Li (2020).
#'
#' When using \code{\link[geostan]{stan_car}}, always use `cmat = TRUE` (the default). 
#'
#' @source
#'
#' Cliff A, Ord J (1981). Spatial Processes: Models and Applications. Pion.
#'
#' Cressie N (2015 \[1993\]). Statistics for Spatial Data. Revised edition. John Wiley & Sons.
#' 
#' Cressie N, Perrin O, Thomas-Agnan C (2005). “Likelihood-based estimation for Gaussian MRFs.” Statistical Methodology, 2(1), 1–16.
#'
#' Cressie N, Wikle CK (2011). Statistics for Spatio-Temporal Data. John Wiley & Sons.
#'
#' Haining RP, Li G (2020). Modelling Spatial and Spatio-Temporal Data: A Bayesian Approach. CRC Press.
#'
#' @return A list containing all of the data elements required by the Stan CAR model.
#'
#' @examples
#'
#' data(georgia)
#'
#' ## binary adjacency matrix
#' A <- shape2mat(georgia, style = "B")
#'
#' ## get list of data for Stan
#' cp <- prep_car_data(A, "WCAR")
#'  # range of permissible values for rho:
#' 1 / range(cp$lambda)
#' 
#' \dontrun{
#' ## pass the data to stan_car
#' fit = stan_car(log(rate.male) ~ 1, data = georgia, car_parts = cp)
#' }
#'
#' # or use the ACAR specification
#' cp <- prep_car_data(A, "ACAR")
#'  # range of permissible values for rho:
#' 1 / range(cp$lambda)
#' @export
#' @md
#' @importFrom rstan extract_sparse_parts
#' @importFrom Matrix isSymmetric Matrix rowSums summary
prep_car_data <- function(A, style = c("WCAR", "ACAR", "DCAR"), k = 1, gamma = 0, lambda = TRUE, cmat = TRUE) {
    style = match.arg(style)
    stopifnot(inherits(A, "matrix") | inherits(A, "Matrix"))
    stopifnot(all(Matrix::rowSums(A) > 0))
    A <- Matrix::Matrix(A, sparse = TRUE)
    n <- nrow(A)
    if (style == "ACAR") {
        Ni <- Matrix::rowSums(A)
        A.idx <- Matrix::summary(A)
        C <- Matrix::Matrix(0, nrow = n, ncol = n)        
        for (m in 1:nrow(A.idx)) {
            i <- A.idx[m, "i"]; j <- A.idx[m, "j"]
            C[i,j] <- sqrt( Ni[j] ) / sqrt( Ni[i] )
        }
        M_diag <- 1 / Ni
    }
    if (style == "DCAR") {
        dinv <- A
        dinv[dinv>0] <- (dinv[dinv>0] + gamma)^(-k)
        max.dinv <- max(dinv)
        dinv <- dinv / max.dinv
        # conditional variance proportional to total d^-k
        dinv.sums <- Matrix::rowSums(dinv)
        M_diag <- 1 / dinv.sums
        # C scaled by sqrt of ratio of total distances
        A.idx <- Matrix::summary(dinv)
        C <- Matrix::Matrix(0, nrow = n, ncol = n)        
        for (m in 1:nrow(dinv.idx)) {
            i <- A.idx[m, "i"]; j <- A.idx[m, "j"]
            C[i,j] <- dinv[i,j] * sqrt(dinv.sums[j] / dinv.sums[i])
        }
    }
    if (style == "WCAR") {
        Ni <- Matrix::rowSums(A)
        C <- A / Ni
        M_diag <- 1 / Ni        
        stopifnot( Matrix::isSymmetric(C %*% Matrix::Diagonal(x = M_diag), check.attributes = FALSE) )
        car.dl <- rstan::extract_sparse_parts(A)
        names(car.dl) <- paste0("Ax_", names(car.dl))
        car.dl$nAx_w <- length(car.dl$Ax_w)
        car.dl$Cidx <- array(0, dim = 1)
        car.dl$nC <- 1
        car.dl$WCAR <- 1
    } else {
        stopifnot( Matrix::isSymmetric(C %*% Matrix::Diagonal(x = M_diag), check.attributes = FALSE) )        
        car.dl <- rstan::extract_sparse_parts(Matrix::Diagonal(n) - C)
        names(car.dl) <- paste0("Ax_", names(car.dl))
        car.dl$nAx_w <- length(car.dl$Ax_w)        
        car.dl$Cidx <- which( car.dl$Ax_w != 1 )
        car.dl$nC <- length(car.dl$Cidx)
        car.dl$WCAR <- 0
    }
    car.dl$Delta_inv <- 1 / M_diag
    car.dl$style <- style
    car.dl$log_det_Delta_inv = base::determinant(diag(car.dl$Delta_inv), log = TRUE)$modulus
    car.dl$n <- n
    if (lambda) {
        MCM <- Matrix::Diagonal(x = 1 / sqrt(M_diag)) %*% C %*% Matrix::Diagonal(x = sqrt(M_diag))
        stopifnot(Matrix::isSymmetric(MCM, check.attributes = FALSE))
        lambda <- eigen(MCM)$values
        cat ("Range of permissible rho values: ", 1 / range(lambda), "\n")
        car.dl$lambda <- lambda
    }
    if (cmat) car.dl$C <- C
    return (car.dl)
}


#' Download shapefiles
#'
#' @description Given a url to a shapefile in a compressed .zip file, download the file and unzip it into a folder in your working directory.
#' 
#' @param url url to download a shapefile.
#' @param folder what to name the new folder in your working directory containing the shapefile
#' 
#' @return A folder in your working directory with the shapefile; filepaths are printed to the console.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' url <- "https://www2.census.gov/geo/tiger/GENZ2019/shp/cb_2019_us_state_20m.zip"
#' folder <- tempdir()
#' print(folder)
#' get_shp(url, folder)
#' states <- sf::st_read(folder)
#' head(states)
#' }
#' @export
#' @importFrom utils unzip download.file
#' 
get_shp <- function(url, folder = "shape") {
	tmp.dir <- tempfile(fileext=".zip")
	download.file(url, destfile = tmp.dir)
	unzip(tmp.dir, exdir = folder)
	list.files(folder, full.names = TRUE)
}

