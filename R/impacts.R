
#' Spillover/diffusion effects for spatial lag models
#' 
#' @export
#' @rdname impacts
#' @md
#'
#' @param beta Coefficient for covariates (numeric vector)
#' 
#' @param gamma Coefficient for spatial lag of covariates 
#' 
#' @param W Spatial weights matrix
#' 
#' @param rho Spatial dependence parameter (single numeric value)
#'
#' @param approx For a computationally efficient approximation to the required matrix inverse (after LeSage and Pace 2009, pp. 114--115); if `FALSE`, then a proper matrix inverse will be computed using `Matrix::solve`.
#'
#' @param K Degree of polynomial in the expansion to use when 'approx = TRUE`.
#'
#' @details
#'
#' These methods apply only to the spatial lag and spatial Durbin lag models (SLM and SDLM) as fit by `geostan::stan_sar`.
#' 
#' The equation for these SAR models specifies simultaneous feedback between all units, such that changing the outcome in one location has a spill-over effect that may extend to all other locations (a ripple or diffusion effect); the induced changes will also react back onto the first unit. (This presumably takes time, even if the observations are cross-sectional.)
#'
#' These spill-overs have to be incorporated into the interpretation and reporting of the regression coefficients of SLM and SDLM models. A unit change in the value of \eqn{X} in one location will impact \eqn{y} in that same place ('direct' impact) and will also impact \eqn{y} elsewhere through the diffusion process ('indirect' impact). The 'total' expected impact of a unit change in \code{X} is the sum of the direct and indirect effects (LeSage and Pace 2009). 
#' 
#' The `spill` function is for quickly calculating average spillover effects given point estimates of parameters. 
#'
#' The `impacts` function calculates the (average) direct, indirect, and total effects once for every MCMC sample to produce samples from the posterior distribution for the impacts; the samples are returned together with a summary of the posterior distribution (mean, median, and select quantiles).
#' 
#' @source
#'
#' LeSage, James and Pace, R. Kelley (2009). *Introduction to Spatial Econometrics*. CRC Press.
#'
#' LeSage, James (2014). What Regional Scientists Need to Know about Spatial Econometrics. *The Review of Regional Science* 44: 13-32 (2014 Southern Regional Science Association Fellows Address).
#'
#' @importFrom Matrix solve rowSums diag
#'
#' @examples
#' ##
#' ## SDLM data
#' ##
#' 
#' parts <- prep_sar_data2(row = 9, col = 9, quiet = TRUE)
#' W <- parts$W
#' x <- sim_sar(w=W, rho=.6)
#' Wx <- (W %*% x)[,1]
#' mu <- .5 * x + .25 * Wx
#' y <- sim_sar(w=W, rho=0.6, mu = mu, type = "SLM")
#' dat <- cbind(y, x)
#'
#' # impacts per the above parameters
#' spill(0.5, 0.25, 0.6, W)
#'
#' ##
#' ## impacts for SDLM
#' ##
#' 
#' fit <- stan_sar(y ~ x, data = dat, sar = parts,
#'                 type = "SDLM", iter = 500,
#'                 slim = TRUE, quiet = TRUE) 
#'
#' # impacts (posterior distribution)
#' impax <- impacts(fit)
#' print(impax)
#' 
#' # plot posterior distributions
#' og = par(mfrow = c(1, 3),
#'          mar = c(3, 3, 1, 1))
#' S <- impax$samples[[1]]
#' hist(S[,1], main = 'Direct')
#' hist(S[,2], main = 'Indirect')
#' hist(S[,3], main = 'Total')
#' par(og)
#'
#' ##
#' ## The approximate method
#' ##
#'
#' # High rho value requires more K; rho^K must be near zero
#' Ks <- c(10, 15, 20, 30, 35, 40)
#' print(cbind(Ks, 0.9^Ks))
#'
#' # understand sensitivity of results to K when rho is high
#' spill(0.5, -0.25, 0.9, W, approx = TRUE, K = 10)
#' spill(0.5, -0.25, 0.9, W, approx = TRUE, K = 20)
#' spill(0.5, -0.25, 0.9, W, approx = TRUE, K = 30)
#' spill(0.5, -0.25, 0.9, W, approx = TRUE, K = 50)
#'
#' # the correct results
#' spill(0.5, -0.25, 0.9, W, approx = FALSE)
#'
#' # moderate and low rho values are fine with smaller K
#' spill(0.5, -0.25, 0.7, W, approx = TRUE, K = 15)
#' spill(0.5, -0.25, 0.7, W, approx = FALSE)
#' 
spill <- function(beta, gamma = 0, rho, W, approx = TRUE, K = 15) {       
    
    if (approx) {
        T <- matrix(0, nrow = 2, ncol = K + 1)
        T[1, 1] <- 1
        for (i in 2:K) T[1, i + 1] <- mean(diag_power(W, i))
        for (i in 2:(K+1)) T[2, i] <- mean(diag_power(W, i))
        return( impacts_multiplier(beta, gamma, rho, T, K) )
    }

    N <- nrow(W)    
    I <- Matrix::diag(N)
    imrw <- I - rho * W
    Multiplier <- Matrix::solve(imrw)                   
    
    M_spills <- Multiplier %*% (I * beta + W * gamma)
    dir = mean(Matrix::diag(M_spills))
    total <- mean(Matrix::rowSums(M_spills))
    indir <- total - dir
    spills = c(direct = dir,
               indirect = indir,
               total = total)    
    return (spills)    
}


#' @param object A fitted spatial lag  model (from `stan_sar`) 
#'
#' @md
#'
#' @export
#' 
impacts <- function(object, approx = TRUE, K = 15) {

    stopifnot(object$spatial$method == "SAR")
    stopifnot(grepl("SLM|SDLM", object$sar_type))
    
    W <- object$C
    rho <- as.matrix(object, "sar_rho")[,1]
    B <- as.matrix(object, "beta")
    Blabs <- colnames(B)
    M <- ncol(B)
    S <- nrow(B)    
    G <- matrix(0, nrow = S, ncol = M)    
    has_gamma <- inherits(object$slx, "formula")
    
    if (has_gamma) {        
        gamma <- as.matrix(object, "gamma")
        gamma_idx <- match( gsub("^w.", "", colnames(gamma)), Blabs )
        for (j in seq_along(gamma_idx)) G[ , gamma_idx[j] ] <- gamma[ , j ]        
    }
    
    impax <- vector("list", length = M)

    if (approx) {
        
        T <- matrix(0, nrow = 2, ncol = K + 1)
        T[1, 1] <- 1
        for (i in 2:K) T[1, i + 1] <- mean(diag_power(W, i))
        for (i in 2:(K+1)) T[2, i] <- mean(diag_power(W, i))

        for (m in 1:M) {
            impax[[m]] <- sapply(1:S, function(s)
                impacts_multiplier(as.numeric( B[s,m] ), as.numeric( G[s,m] ), rho[s], T, K)) |>
                t()            
        }
        
    } else {

        for (m in 1:M) {
            impax[[m]] <- sapply(1:S, function(s)
                spill(beta = as.numeric( B[s, m] ),
                      gamma = as.numeric( G[s, m] ),
                      rho = rho[s],
                      W = W,
                      approx = approx,
                      K = K)
                ) |>
                t()                  
        }                
    }
    
    summary <- vector("list", length = M)
    for ( m in 1:M) {
        est = apply(impax[[m]], 2, mean)
        est2 = apply(impax[[m]], 2, median)
        sd = apply(impax[[m]], 2, sd)
        lwr = as.numeric(apply(impax[[m]], 2, quantile, probs = 0.025))
        upr = as.numeric(apply(impax[[m]], 2, quantile, probs = 0.975))
        res <- cbind(mean = est, median = est2, sd, lwr, upr)
        row.names(res) <- c('Direct', 'Indirect', 'Total')
        summary[[m]] <- res 
    }

    names(impax) <- Blabs
    names(summary) <- Blabs
    matrix_inverse_method <- ifelse(approx, 'approximate', 'proper')
    out <- list(summary = summary, samples = impax, method = matrix_inverse_method)
    class(out) <- append("impacts_slm", class(out))    
    return(out)    
}

#' @export
#' 
#' @param x An object of class 'impacts_slm', as returned by `geostan::impacts` 
#'
#' @param digits Round results to this many digits
#'
#' @param ... Additional arguments will be passed to `base::print`
#' 
#' @rdname impacts
print.impacts_slm <- function(x, digits = 2, ...) {
    print(x$summary, digits = digits, ...)
}


#' After LeSage and Pace 2009 pp. 114--115
#' @noRd
impacts_multiplier <- function(beta, gamma, rho, T, K) {
                     
    g = numeric(K + 1)
    for (i in 0:K) g[i + 1] <- rho^i
    G <- diag(g)
    
    P <- cbind(beta, gamma)
    
    # direct
    direct = sum(P %*% T %*% G)
    
    # total
    total <- (beta + gamma) * sum(g)
    
    # indirect
    indirect <- total - direct

    return (c(direct = direct, indirect = indirect, total = total))
}

#' diagonal entries of matrix powers e..g, diag( W^{20} )
#' @noRd
diag_power <- function(W, K = 20, trace = FALSE) {
    N <- nrow(W)    
    ones <- matrix(1, nrow = N, ncol = 1)
    P <- K - 1
    mat <- W
    if (K > 2) {
        for (j in 2:P) {
            mat <- W %*% mat
        }
    }
    dwk <- as.numeric( (W * t(mat)) %*% ones )
    if (trace) return (sum(dwk))
    return (dwk)
}

