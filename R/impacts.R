

spill <- function(b1, b2 = 0, W, rho, short = FALSE, K = 4) {
    
    N <- nrow(W)
    I <- diag(rep(1, N))

    if (short == TRUE) {
        
        warning('this does not return accurate results.')

        M_tmp <- I + rho * W
        W_power <- W    
        k <- 2
        while (k >= K) {
            rk <- rho^k
            W_power <- W_power %*% W
            M_tmp <- M_tmp + rk * W_power
            k <- K + 1
        }
        
    } else {
        
        imrw <- I - rho * W
        M_tmp <- Matrix::solve(imrw)
        
    }
    
    M <- M_tmp %*% (I * b1 + W * b2)
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
#' @param short Used short-cut matrix inverse? Not recommended currently.
#'
#' @param K For the short-cut method. Number of powers to use in the Taylor expansion for the short-cut matrix inverse.
#'
#' @param samples Return MCMC samples together with summary in a list? If \code{FALSE}, only the summary is returned.
#'
#'
#' @md
#'
#' @export
#' 
impacts <- function(object, W = object$C, short = FALSE, K = 4, samples = TRUE) {
    
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
                    B[i,m],
                    G[i,Gidx[m]],
                    W = W,
                    rho = rho[i],
                    short = short,
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


