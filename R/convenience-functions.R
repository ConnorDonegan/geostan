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
aple <- function(x, w, digits = 3) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix"))
    stopifnot(all(dim(w) == length(x)))    
    if (any(rowSums(w) != 1)) {
        message("Row standardizing w with: w <- w / rowSums(w)")
        w <- w / rowSums(w)
    }    
    z <- as.numeric(scale(x))
    n <- length(z)
    I <- diag(n)
    lambda <- eigen(w)$values
    w2 <- (w + t(w)) / 2
    top <- t(z) %*% w2 %*% z
    wl <- (t(w) %*% w + as.numeric(t(lambda) %*% lambda) * I / n)
    bottom <- t(z) %*% wl %*% z
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
#' @param rho Spatial autocorrelation parameter in the range [-1, 1]. Typically a scalar value; otherwise an n-length numeric vector.
#' @param sigma Scale parameter (standard deviation). Defaults to \code{sigma = 1}. Typically a scalar value; otherwise an n-length numeric vector.
#' @param ... further arguments passed to \code{MASS::mvrnorm}.
#' 
#' @return
#'
#' If \code{m = 1} a vector of the same length as \code{mu}, otherwise an \code{m x length(mu)} matrix with one sample in each row.
#'
#' @details Calls \code{MASS::mvrnorm} internally to draw from the multivariate normal distribution. The covariance matrix is specified following the simultaneous autoregressive (SAR) model. 
#'
#' @seealso \link[geostan]{aple}, \link[geostan]{shape2mat}, \link[MASS]{mvrnorm} 
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
    stopifnot(inherits(w, "matrix"))
    stopifnot(ncol(w) == nrow(w))
    stopifnot(all(rowSums(w) %in% c(0, 1)))
    N <- nrow(w)
    if (missing(mu)) {
        mu <- rep(0, N)
    } else {
        if (length(mu) != N | !inherits(mu, "numeric")) stop("mu must be a numeric vector with length equal to nrow(W).")
        }
    if (!inherits(rho, "numeric") | !length(rho) %in% c(1, N) | any(rho > 1 | rho < -1)) stop("rho must be numeric value within range [-1, 1], or a k-length numeric vector where K=nrow(W).")
    if (!inherits(sigma, "numeric") | !length(sigma) %in% c(1, N) | !all(sigma > 0)) stop("sigma must be a positive numeric value, or k-length numeric vector, with K=nrow(W).")
    I <- diag(N)
    S <- sigma^2 * solve( (I - rho * t(w)) %*% (I - rho * w) )
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
#' @details If any observations with no neighbors are found (i.e. \code{any(rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
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
#'
#' 
mc <- function(x, w, digits = 3, warn = TRUE) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix"))
    stopifnot(all(dim(w) == length(x)))
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
        if (warn) message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data.")
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
#' If any observations with no neighbors are found (i.e. \code{any(rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
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
#' library(sf)
#' data(georgia)
#' x <- georgia$ICE
#' w <- shape2mat(georgia, "W")
#' moran_plot(x, w)
#'
#' @import ggplot2
#' @importFrom signs signs
moran_plot <- function(x, w, xlab = "x (centered)", ylab = "Spatial Lag", pch = 20, col = "darkred", size = 2, alpha = 1, lwd = 0.5) {
    stopifnot(inherits(x, "numeric") | inherits(x, "integer"))
    stopifnot(inherits(w, "matrix"))
    stopifnot(all(dim(w) == length(x)))
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
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
#' @param w An `n x n` spatial connectivity matrix. See \link[geostan]{shape2mat}. If \code{w} is not row standardized (\code{all(rowSums(w) == 1)}), it will automatically be row-standardized.
#' @param type Return the type of association also (High-High, Low-Low, High-Low, and Low-High)? Defaults to \code{FALSE}.
#'
#' @details
#'
#' The values will be standardized with \code{z = scale(x)} first and \code{w} will be row-standardized if needed. The LISA values are the product of each \code{z} value with their respective mean surrounding value \code{lagz = w \%*\% z}; \code{lisa = z * lagz}. These are for exploratory analysis and model diagnostics. The function uses Equation 7 from Anselin (1995).
#'
#' An above-average value (i.e. positive z-value) with positive mean spatial lag indicates local positive spatial autocorrelation and is designated type "High-High"; a low value surrounded by high values indicates negative spatial autocorrelation and is designated type "Low-High", and so on.
#' 
#' @return If \code{type = FALSE} a numeric vector of lisa values for exploratory analysis of local spatial autocorrelation. If \code{type = TRUE}, a \code{data.frame} with columns \code{zi} (the lisa value) and \code{type}.
#'
#' @seealso \link[geostan]{moran_plot}, \link[geostan]{mc}, \link[geostan]{aple}
#'
#' @source
#'
#' Anselin, Luc. "Local indicators of spatial association—LISA." Geographical analysis 27, no. 2 (1995): 93-115.
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
#' ggplot(georgia, aes(fill = li)) +
#'   geom_sf() +
#'   scale_fill_gradient2()
#' 
lisa <- function(x, w, type = FALSE) {
    stopifnot(length(x) == nrow(w) & length(x) == ncol(w))
    if (any(rowSums(w) != 1)) {
        message("Row standardizing w with: w <- w / rowSums(w)")
        w <- w / rowSums(w)
    }    
    z <- scale(x)
    lag <- as.numeric(w %*% z)
    zi <- as.numeric(z * lag)
    if (any(rowSums(w) == 0)) {
        zero.idx <- which(rowSums(w) == 0)
        zi[zero.idx] <- NA        
        message(length(zero.idx), " observations with no neighbors found. They will be converted to NA.")
    }     
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
#' @seealso \link[geostan]{me_diag}, \link[geostan]{mc}, \link[geostan]{moran_plot}, \link[geostan]{aple}
#' 
#' @export
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(sf)
#' data(georgia)
#' sp_diag(georgia$college, georgia)
#' fit <- stan_glm(college ~ 1, data = georgia, refresh = 0)
#' sp_diag(fit, georgia)
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
    stopifnot(length(y) == nrow(shape))
    shape$y <- y
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
        return( gridExtra::grid.arrange(hist, g.mc, map.y, ncol = 3) )
    } else return (list(hist, mc.y, map.y))
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
#' @param index Integer value; use this if you wish to identify observations with the largest `n=index` absolute Delta values; data on the top `n=index` observations ordered by absolute Delta value will be printed to the console and the plots will be labeled with the indices of the identified observations.
#' 
#' @return A grid of spatial diagnostic plots for measurement error (i.e. data) models comparing the raw observations to the posterior distribution of the true values (given the specified observations, standard errors, and model). The Moran scatter plot and map depict the difference between the posterior means and the raw observations (i.e. shrinkage). 
#'
#' @seealso \link[geostan]{sp_diag}, \link[geostan]{moran_plot}, \link[geostan]{mc}, \link[geostan]{aple}
#' @export
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
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
#' ## see index values for the largest delta values
#'  ## (differences between raw estimate and the posterior mean)
#' me_diag(fit, "ICE", georgia, index = 3)
#' }
#' 
me_diag <- function(fit,
                    varname,                    
                    shape,
                    w,
                    probs = c(0.025, 0.975),
                    plot = TRUE,
                    size = 0.25,
                    index = 0
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
    df$Delta <- x.mu - x.raw
    if (plot) {
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
            labs(
                x = xlbl,
                y = paste("Posterior Mean and", width, "C.I.")
            ) +                
            theme_classic()
        if (!missing(w)) {
            g.mc <- moran_plot(df$Delta, w) +
                labs(
                    x = expression(paste(Delta, ' (centered)'))
                     )
        } else {
            w <- shape2mat(shape, style = "W")
            g.mc <- moran_plot(df$Delta, w) +
                labs(x = expression(paste(Delta, ' (centered)')))            
        }
        map.delta <- ggplot(shape) +
            geom_sf(aes(fill = df$Delta),
                    lwd = 0.05,
                    col = "gray20") +
            scale_fill_gradient2(name = bquote(Delta)) +
            theme_void() +
            theme(
                legend.position = "bottom"
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
            g.mc <- g.mc +
                geom_label(aes(label = df$ID),
                           alpha=0.9,
                           na.rm = TRUE
                           )
            map.delta <- map.delta +
                geom_sf(aes(col = !is.na(df$ID),
                            fill = NULL),
                        fill = alpha("white", 0)) +
                scale_color_manual(
                    breaks =  c(FALSE, TRUE),
                    values = c(NA, "black"),
                    guide = FALSE
                    )
                    
        }
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
#' @description A wrapper function for a string of \link{spdep} (and other) functions required to convert spatial objects to spatial or spatio-temporal connectivity matrices.
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
#' @seealso \code{\link[geostan]{edges}}
#' 
#' @details
#'
#' Haining and Li (Ch. 4) provide a helpful discussion of spatial connectivity matrices (Ch. 4) and spatio-temporal connectivity matricies (Ch. 15).
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
#' 
#' data(georgia)
#'
#' ## binary adjacency matrix
#' C <- shape2mat(georgia, "B")
#' ## row sums gives the numbers of neighbors per observation
#' rowSums(C)
#' diag(C)
#'
#' ## row-standardized matrix 
#' W <- shape2mat(georgia, "W")
#' rowSums(W)
#' diag(W)
#' 
#' ## for space-time data
#' ## if you have multiple years with same neighbors
#' ## provide the geography (for a single year!) and number of years \code{t}
#' ## defaults to the contemporaneous connectivity structure
#' Cst <- shape2mat(georgia, t = 5)
#' 
shape2mat <- function(shape, style = c("B", "W", "C", "S"), t = 1, st.type = "contemp", zero.policy = TRUE, queen = TRUE, snap = sqrt(.Machine$double.eps)) {
  style <- match.arg(style)
  shape_class <- class(shape)
  if (!any(c("sf", "SpatialPolygonsDataFrame", "SpatialPolygons") %in% shape_class)) stop("Shape must be of class SpatialPolygons, SpatialPolygonsDataFrame, or sf (simple features).")
      w <- spdep::nb2mat(spdep::poly2nb(shape, queen = queen, snap = snap), style = style, zero.policy = zero.policy)
  attributes(w)$dimnames <- NULL
  if (t > 1) { 
      if (style != "B") stop ("Only the binary coding scheme (style = 'B') has been implemented for space-time matrices, but you can row-standardize any matrix using: 'W = C / rowSums(C)' where C is a binary-coded connectivity matrix.")
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
#' @description Widely Application Information Criteria (WAIC) for model evaluation
#' @export
#' @param fit An \code{geostan_fit} object or any Stan model with a parameter named "log_lik", the pointwise log predictive likelihood
#' @param pointwise Logical, should a vector of values for each observation be returned? 
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
#'  Chun, Yongwan and Griffith, Daniel A. (2013). Spatial statistics and geostatistics. Sage, p. 18.
#' 
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
#' 
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. Spatial Statistics. \url{https://doi.org/10.1016/j.spasta.2020.100450}
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
#' 
#' @return
#' 
#' Returns a \code{data.frame} with three columns. The first two columns (\code{node1} and \code{node2}) contain the indices of connected pairs of nodes; only unique pairs of nodes are included. The third column (\code{weight}) contains the element \code{w[node1, node2]}.
#'
#' @details This is used internally for \link[geostan]{stan_icar} and it is also helpful for creating the scaling factor for BYM2 models fit with \code{stan_icar}.
#'
#' @seealso \link[geostan]{shape2mat}, \link[geostan]{prep_icar_data}, \link[geostan]{stan_icar}
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
    if(length(nbs)) data.frame(node1 = i, node2 = nbs, weight = w[i, nbs])
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
#' @description Given a symmetric n x n connectivity matrix, prepare data for intrinsic conditional autoregressive models in Stan.
#' 
#' @param C connectivity matrix
#' @param inv_sqrt_scale_factor optional vector of scale factors for each connected portion of the graph structure. If not provided by the user it will be fixed to a vector of ones. 
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
#' \item{inv_sqrt_scale_factor}{k-length vector of ones. Placeholder for user-specified information.}
#' \item{comp_id}{n-length vector indicating the group membership of each observation.}
#' }
#'
#' @details
#'
#' This is used internally to prepare data for \link[geostan]{stan_icar} models. It can also be helpful for fitting custom ICAR models outside of \code{geostan}. It relies on \link[spdep]{graph2nb} and \link[spdep]{n.comp.nb}.
#' 
#' @seealso \link[geostan]{stan_icar}, \link[geostan]{edges}, \link[geostan]{shape2mat}
#'
#' @examples
#' 
#' data(sentencing)
#' C <- shape2mat(sentencing)
#' icar.data.list <- prep_icar_data(C)
#' 
#' @export
#' @importFrom spdep n.comp.nb graph2nb
prep_icar_data <- function (C, inv_sqrt_scale_factor = NULL) {
  n <- nrow(C)
  E <- edges(C)
  G <- list(np = nrow(C), from = E$node1, to = E$node2, nedges = nrow(E))
  class(G) <- "Graph"
  nb2 <- spdep::n.comp.nb(spdep::graph2nb(G))
  k = nb2$nc
  if (inherits(inv_sqrt_scale_factor, "NULL")) inv_sqrt_scale_factor <- array(rep(1, k), dim = k)
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
    A <- model.matrix(~ 0, data.frame(C))
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
            inv_sqrt_scale_factor = inv_sqrt_scale_factor, 
            comp_id = nb2$comp.id)
  return(l)
}


#' Prepare data for Stan CAR model
#'
#' @export
#' @importFrom rstan extract_sparse_parts
#' 
#' @param A Binary adjacency matrix.
#' 
#' @param style Specification for the connectivity matrix (C) and conditional variances (M); one of "WCAR", "ACAR", or "DCAR".
#' 
#' @param lambda If TRUE, return eigenvalues required for calculating the log determinant of the precision matrix and for determining the range of permissible values of rho. These will also be printed with a message if lambda = TRUE.
#' 
#' @param cmat Return the full matrix C if TRUE.
#' 
#' @param d distance matrix. Non-neighboring values (as indicated by A) will be set to zero internally.
#' 
#' @param k If provided with `style = DCAR` and `d`, distances will be raised to the -k power (d^-k).
#' 
#' @details The CAR model is N(Mu, Sigma), Sigma = (I - rho C)^-1 M.
#'
#' The DCAR specification is distance-based, and requires the user provide a distance matrix \code{d} in addition to adjacency matrix \code{A}. The WCAR specification uses row-standardized connectivity matrix and sets conditional variances proportional to the number of neighbors. The ACAR specification is from Cressie, Perrin and Thomas-Agnon (2005); also see Cressie and Wikle (2011, p. 188).
#'
#' For inverse-distance weighting schemes, see Cliff and Ord (1981); for distance-based CAR specifications, see Cressie (2015 [1993]) and Haining and Li (2020).
#'
#' \code{stan_car} requires the user to provide the full C matrix (which is returned when cmat = TRUE).
#'
#' @source
#'
#' Cliff A, Ord J (1981). Spatial Processes: Models and Applications. Pion.
#'
#' Cressie N (2015 [1993]). Statistics for Spatial Data. Revised edition. John Wiley & Sons.
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
#'
#' ## pass the data to stan_car
#' fit = stan_car(ICE ~ 1, data = georgia, car_parts = cp)
#'
#' ## diagnostics
#' sp_diag(fit, georgia)
#' 
prep_car_data <- function(A, style = c("ACAR", "WCAR", "DCAR"), lambda = FALSE, cmat = TRUE, d, k = 1) {
    style = match.arg(style)
    n <- nrow(A)    
    if (style == "WCAR") {
        Ni <- rowSums(A)
        C <- A / Ni
        M_diag <- 1 / Ni
    }
    if (style == "ACAR") {
        Ni <- rowSums(A)
        C <- matrix(0, nrow = n, ncol = n)
        for (i in 1:n) for (j in 1:n) C[i,j] <- A[i,j ] * sqrt( Ni[j] ) / sqrt( Ni[i] )
        M_diag <- 1 / Ni
    }
    if (style == "DCAR") {
        if (missing(d)) stop("d matrix is required for the DCAR specification.")
        class(d) <- "matrix"    
        stopifnot(all( dim(d) == n ))        
        # d -> d^-k; neighbors only
        dinv <- A * d
        dinv[dinv>0] <- dinv[dinv>0]^(-k)
        # scale by max d^-k
        max.dinv <- max(dinv)
        dinv <- dinv / max.dinv
        # conditional variance proportional to total d^-k
        dinv.sums <- rowSums(dinv)
        M_diag <- 1 / dinv.sums
        # C scaled by ratio of conditional standard deviations
        C <- matrix(0, nrow = n, ncol = n)
        for (i in 1:n) for (j in 1:n) C[i,j] <- dinv[i,j] * sqrt( dinv.sums[j] / dinv.sums[i]  )
    }
    stopifnot( isSymmetric.matrix(C %*% diag(M_diag), check.attributes = FALSE) )
    car.dl <- rstan::extract_sparse_parts(diag(n) - C)
    car.dl$Cidx <- which( car.dl$w != 1 )
    names(car.dl)[1] <- "ImC"
    names(car.dl)[2] <- "ImC_v"
    names(car.dl)[3] <- "ImC_u"
    car.dl$nImC <- length(car.dl$ImC)
    car.dl$nC <- length(car.dl$Cidx)
    car.dl$M_diag <- M_diag
    car.dl$dim_C <- n
    car.dl$style <- style
    if (lambda) {
        MCM <- diag( 1 / sqrt(M_diag) ) %*% C %*% diag( sqrt(M_diag) )
        stopifnot(isSymmetric.matrix(MCM, check.attributes = FALSE))
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

#' Theil's inequality index
#'
#' @export
#' 
#' @md
#' 
#' @description Calculates Theil's entropy-based measure of inequality for a given set of disease incidence rates and populations at risk. 
#'
#' @param count Number of cases (e.g., disease incidence). 
#' @param population Population size (population at risk).
#' @param rate Incidence rates. If provided, case counts will be calculated automatically as `cases = rates * Population`#' 
#' @param total If `TRUE`, the values will all be summed; if `FALSE`, then each area's contribution to total inequality will be returned.
#' @return if `total = TRUE`, a scalar value; if `total = FALSE`, a vector of numeric values, where each value represents that area's contribution to total inequality.
#' 
#' @details
#'
#' Theil's index is a good index of inequality in disease and mortality burdens when multiple groups are being considered, as is typical of geospatial analysis. It provides a summary measure of inequality across a set of areal units, such as counties, that may be tracked over time. Also, it is interesting because it is additive, and thus admits of simple decompositions. 
#'
#' The index measures discrepancies between a population's share of the disease burden, `omega`, and their share of the population, `eta`. A situation of zero inequality would imply that each population's share of cases is equal to its population share, or, `omega=eta`. Each population's contribution to total inequality is calculated as:
#' ```
#'              T_i = omega_i * [log(omega_i/eta_i)],
#' ```
#' the log-ratio of case-share to population-share, weighted by their share of cases. Theil's index for all areas is the sum of each area's T_i:
#' ```
#'              T = sum_(i=1)^n T_i.
#' ```
#' Theil's T is thus a weighted mean of log-ratios of case shares to population shares, where each log-ratio (which we may describe as a raw inequality score) is weighted by its share of total cases. The index has a minimum of zero and a maximum of `log(N)`, where `N` is the number of units (e.g., number of counties).

#' Theil's index is based on Shannon's information theory, he used it to study a variety of topics, including income inequality and racial segregation. Theil's index is often of great interest because it is additive across multiple scales, such as when the data has a nested structure to it (e.g., counties within states). The Texas Inequality Project provides introductions to, and examples of using, the Theil index (Conceicao and Ferreira, 2000). However, this `R` function is just a simple implementation for `flat' or non-nested data structures (e.g., a set of counties). 
#'
#' @source
#'
#' Conceicao, P. and P. Ferreira (2000). The young person's guide to the Theil Index: Suggesting intuitive interpretations and exploring analytical applications. University of Texas Inequality Project. UTIP Working Paper Number 14. Accessed May 1, 2021 from \url{https://utip.gov.utexas.edu/papers.html}
#'
#' Theil, Henri (1972). *Statistical Decomposition Analysis.* Amsterdan, The Netherlands and London, UK: North-Holland Publishing Company.
#'
#' Shannon, Claude E. and Weaver, Warren (1963). *The Mathematical Theory of Communication*. Urbana and Chicago, USA: University if Illinois Press.
#'
#' @examples
#'
#'
theil <- function(count, population, rates, total = TRUE) {
    if (missing(count)) count <- rates * population
    omega = count / sum(count)
    eta = population / sum( population )
    T = omega * log (omega / eta)
    T[is.na(T)] <- 0
    if (total) T = sum( T )
    return (Theil = T)
}

