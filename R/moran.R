
#' The Moran coefficient
#'
#' @description The Moran coefficient, a measure of spatial autocorrelation (also known as Global Moran's I)
#' 
#' @export
#' 
#' @param x Numeric vector of input values, length n.
#' @param w An n x n spatial connectivity matrix. See \code{\link[geostan]{shape2mat}}. 
#' @param digits Number of digits to round results to.
#' @param warn If `FALSE`, no warning will be printed to inform you when observations with zero neighbors or `NA` values have been dropped. 
#' @param na.rm If `na.rm = TRUE`, observations with `NA` values will be dropped from both `x` and `w`. 
#' @return The Moran coefficient, a numeric value.
#'
#' @details If any observations with no neighbors are found (i.e. \code{any(Matrix::rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#'
#' @seealso \link[geostan]{moran_plot}, \link[geostan]{lisa}, \link[geostan]{aple}
#' 
#' @examples
#' library(sf)
#' data(georgia)
#' w <- shape2mat(georgia, style = "W")
#' x <- georgia$ICE
#' mc(x, w)
#' @source
#'
#' Chun, Yongwan, and Daniel A. Griffith. Spatial Statistics and Geostatistics: Theory and Applications for Geographic Information Science and Technology. Sage, 2013.
#' 
#' Cliff, Andrew David, and J. Keith Ord. Spatial processes: models & applications. Taylor & Francis, 1981.
#' 
#' @importFrom Matrix rowSums
#' 
mc <- function(x, w, digits = 3, warn = TRUE, na.rm = FALSE) {
    check_sa_data(x, w) 
    na_idx <- which(is.na(x))   
    if (na.rm == TRUE && length(na_idx) > 0) {   
        if (warn) message(length(na_idx), " NA values found in x. They will be removed from the data before calculating the Moran coefficient. If matrix w was row-standardized, it may not longer be. To address this, you can use a binary connectivity matrix, using style = 'B' in shape2mat.")
        x <- x[-na_idx]
        w <- w[-na_idx, -na_idx]
    }
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        if (warn) message(length(zero.idx), " observations with no neighbors found. They will be removed from the data before calculating the Moran coefficient.")
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
#' @param na.rm If `na.rm = TRUE`, any observations of `x` with `NA` values will be dropped from `x` and from `w`.
#' @details For details on the symbol parameters see the documentation for \link[ggplot2]{geom_point}.
#'
#' If any observations with no neighbors are found (i.e. \code{any(Matrix::rowSums(w) == 0)}) they will be dropped automatically and a message will print stating how many were dropped.
#' 
#' @return Returns a \code{gg} plot, a scatter plot with \code{x} on the horizontal and its spatially lagged values on the vertical axis (i.e. a Moran scatter plot).
#'
#' @seealso \link[geostan]{mc}, \link[geostan]{lisa}, \link[geostan]{aple}
#'
#' @source
#'
#' Anselin, Luc. "Local indicators of spatial association—LISA." Geographical analysis 27, no. 2 (1995): 93-115.
#' 
#' @examples
#' data(georgia)
#' x <- georgia$income
#' w <- shape2mat(georgia, "W")
#' moran_plot(x, w)
#' @import ggplot2
#' @importFrom signs signs
#' @importFrom Matrix rowSums
moran_plot <- function(x, w, xlab = "x (centered)", ylab = "Spatial Lag", pch = 20, col = "darkred", size = 2, alpha = 1, lwd = 0.5, na.rm = FALSE) {
    check_sa_data(x, w)
    na_idx <- which(is.na(x))
    if (length(na_idx) > 0) {
        if (na.rm == TRUE) {   
            message(length(na_idx), " NA values found in x. They will be dropped from the data before creating the Moran plot. If matrix w was row-standardized, it no longer is. To address this, you can use a binary connectivity matrix, using style = 'B' in shape2mat.")
            x <- x[-na_idx]
            w <- w[-na_idx, -na_idx]
        } else {
            stop("NA values found in moran_plot(x, ...). Use moran_plot(x, na.rm = TRUE, ...)")
        }
    }
    if (any(Matrix::rowSums(w) == 0)) {
        zero.idx <- which(Matrix::rowSums(w) == 0)
        message(length(zero.idx), " observations with no neighbors found. They will be dropped from the data before creating the Moran plot. If matrix w was row-standardized, it no longer is. To address this, you can use a binary connectivity matrix, using style = 'B' in shape2mat.")
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
#' @description A local indicator of spatial association (LISA) based on Moran's I (the Moran coefficient) for exploratory data analysis. The interpretation of LISA values depends on the type of association, which may be negative or positive. However, local Moran's I is best suited for detecting clustering behavior (positive autocorrelation).
#'
#' @param x Numeric vector of length `n`.
#' @param w An `n x n` spatial connectivity matrix. See \link[geostan]{shape2mat}. If \code{w} is not row standardized (\code{all(Matrix::rowSums(w) == 1)}), it will automatically be row-standardized.
#' @param type Return the type of association also (High-High, Low-Low, High-Low, and Low-High)? Defaults to \code{FALSE}.
#'
#' @details
#'
#' The values of `x` will be standardized with \code{z = scale(x)} first and \code{w} will be row-standardized if needed. The LISA values are the product of each \code{z} value with their respective mean surrounding value \code{lagz = w \%*\% z}; \code{lisa = z * lagz}. These are for exploratory analysis and model diagnostics. The function uses Equation 7 from Anselin (1995).
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
#' library(ggplot2)
#' library(sf)
#' data(georgia)
#' w <- shape2mat(georgia, "W")
#' x <- georgia$ICE
#' li = lisa(x, w)
#' head(li)
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
