
#' @export
#' @rdname impacts
#' @md
#' @importFrom Matrix solve rowSums diag
#'
#' @param beta Coefficient for covariates (numeric vector)
#' @param gamma Coefficient for spatial lag of covariates (numeric vector)
#' @param W Spatial weights matrix
#' @param rho Estimate of spatial dependence parameter
#' @param quick Use short-cut method to calculate inverse of matrix (I - rho * W)?
#' @param K Degree of polynomial in the expansion to use for the 'quick' matrix inverse method.
#'
#' @details
#' The `spill` function is for quickly calculating spillover effects given point estimates of parameters. This is used internally by `impacts`, and can be applied iteratively to a matrix of MCMC samples (for beta, gamma, and rho) to obtain MCMC samples for impacts.
#' 
spill <- function(beta, gamma = 0, W, rho, quick = FALSE, K = 15) {
    
    N <- nrow(W)
    I <- diag(rep(1, N))

    if (short == TRUE) {        

        M_tmp <- I + rho * W
        W_k <- W    
        for (j in 2:K) {
            W_k <- W %*% W_k
            M_tmp = M_tmp + rho^j * W_k
        }        
        
    } else {
        
        imrw <- I - rho * W
        M_tmp <- Matrix::solve(imrw)
        
    }
    
    M <- M_tmp %*% (I * beta + W * gamma)
    dir = mean(Matrix::diag(M))
    total <- mean(Matrix::rowSums(M))
    indir <- total - dir
    spills = c(dir, indir, total)    
    return (spills)
    
}

#' Calculate local and global spillover effects
#'
#' @description
#'
#' This needs updating to handle SLM which has 'impacts'
#' but has no slx terms.
#'
#' @param object A fitted geostan SAR model that has spillover effects
#'
#' @param W Spatial weights matrix used to fit the model
#'
#' @param quick Used short-cut matrix inverse? 
#'
#' @param K For the 'quick' method. Number of powers to use in the Taylor expansion for the short-cut matrix inverse.
#'
#' @param samples Return MCMC samples together with summary in a list? If \code{FALSE}, only the summary is returned.
#'
#' @md
#'
#' @export
#' 
impacts <- function(object, W = object$C, quick = FALSE, K = 6, samples = TRUE) {
    
    stopifnot(object$spatial$method == "SAR")
    B <- as.matrix(object, "beta")
    G <- as.matrix(object, "gamma")
    rho <- as.matrix(object, "sar_rho")[,1]
    Blabs <- colnames(B)
    Gidx <- match( gsub("^w.", "", colnames(G)), Blabs )
    M <- length(Gidx) 
    S <- nrow(B)    
    impax <- vector("list", length = M)

    for (m in 1:M) {
        for (i in 1:S) {
            impax[[m]] <- rbind(
                impax[[m]],
                spill(
                    beta = B[i,m],
                    gamma = G[i,Gidx[m]],
                    W = W,
                    rho = rho[i],
                    quick = quick,
                    K = K)
            )               
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
        names(summary)[m] <- Blabs[Gidx[m]]
    }
    
    if (samples == TRUE) return(list(summary = summary, samples = impax))
    
    return(summary)
    
}


