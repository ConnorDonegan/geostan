
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
#' @param method Use 'quick' for a computationally efficient approximation (after LeSage and Pace 2009, pp. 114--115) and 'proper' to compute the matrix inversion using `Matrix::solve`.
#'
#' @param K Degree of polynomial in the expansion to use for the 'quick' method.
#'
#' @details
#'
#' These functions are for interpreting the results of the spatial lag model ('SLM') and the spatial Durbin lag model ('SDLM'). The models can be fitted using \link[geostan]{stan_sar}. The equation for these SAR models specifies simultaneous feedback between all units, such that changing the outcome in one location has a spill-over effect that may extend to all other locations (a ripple or diffusion effect); the induced changes will also react back onto the first unit. (This presumably takes time, even if the observations are cross-sectional.)
#'
#' These spill-overs have to be incorporated into the interpretation of the regression coefficients of SLM and SDLM models (granting that the model specification itself is reasonable for your application). A unit change in the value of \code{X} in one location will impact \code{y} in that same place ('direct' impact) and will also impact \code{y} elsewhere through the diffusion process ('indirect' impact). The 'total' impact of a unit change in \code{X} is the sum of the direct and indirect effects (LeSage and Pace 2009). 
#' 
#' The `spill` function is for quickly calculating average spillover effects given point estimates of parameters. 
#'
#' The `impacts` function calculates the (average) direct, indirect, and total effects once for every MCMC sample to produce samples from the posterior distribution for the impacts; the samples are returned together with a summary of the posterior distribution (mean, median, and select quantiles).
#'
#' These methods are only required for the spatial lag and spatial Durbin lag models (SLM and SDLM).
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
#' impax$summary
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
#' 
spill <- function(beta, gamma = 0, rho, W, method = c('quick', 'proper'), K = 20) {
    
    method <- match.arg(method)
    
    if (method == 'quick') {
        T <- matrix(0, nrow = 2, ncol = K + 1)
        T[1, 1] <- 1
        for (i in 2:K) T[1, i + 1] <- mean(diag_power(W, i))
        for (i in 2:(K+1)) T[2, i] <- mean(diag_power(W, i))
        return( impacts_multiplier(beta, gamma, rho, T, K) )
    }

    # matrix powers: slower method of approximation
    ## M_tmp <- I + rho * W
    ## W_k <- W
    ## for (j in 2:K) {
    ##     W_k <- W %*% W_k
    ##     M_tmp = M_tmp + rho^j * W_k

    N <- nrow(W)    
    I <- diag(rep(1, N))
    imrw <- I - rho * W
    M_tmp <- Matrix::solve(imrw)                   
    
    M <- M_tmp %*% (I * beta + W * gamma)
    dir = mean(Matrix::diag(M))
    total <- mean(Matrix::rowSums(M))
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
impacts <- function(object, method = c('quick', 'proper'), K = 20) {

    stopifnot(object$spatial$method == "SAR")
    stopifnot(grepl("SLM|SDLM", object$sar_type))
    method <- match.arg(method)
    
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

    if (method == "proper") {
        
        for (m in 1:M) {
                impax[[m]] <- sapply(1:S, function(s)
                    spill(beta = B[s, m],
                          gamma = G[s, m],
                          rho = rho[s],
                          W = W,
                          method,
                          K = K)
                    ) |>
                    t()                  
        }
        
    } else {
        
        T <- matrix(0, nrow = 2, ncol = K + 1)
        T[1, 1] <- 1
        for (i in 2:K) T[1, i + 1] <- mean(diag_power(W, i))
        for (i in 2:(K+1)) T[2, i] <- mean(diag_power(W, i))

        for (m in 1:M) {
            impax[[m]] <- sapply(1:S, function(s)
                impacts_multiplier(B[s,m], G[s,m], rho[s], T, K)) |>
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
        names(summary)[m] <- Blabs[m]        
    }    
    
    return(list(summary = summary, samples = impax))
    
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
