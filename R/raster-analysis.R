#' Eigenvalues of a spatial weights matrix: for spatial regression with raster data
#' 
#' @description Approximate eigenvalues for the row-standardized spatial connectivity matrix W of a regular tessellation, e.g., remotely sensed imagery.
#' 
#' @param row Number of rows in the raster dataset
#' @param col Number of columns in the raster dataset
#'
#' @details
#'
#' Uses Equation 5 from Griffith (2000) to calculate the eigenvalues for a row-standardized spatial weights matrix; this is valid for a regular tessellation (rectangular grid or raster). The rook criteria is used to define adjacency.
#'
#' The purpose is to calculate eigenvalues of the spatial weights matrix for the CAR and SAR models, enabling spatial regression with large raster data sets. This function is used internally by \code{\link[geostan]{prep_sar_data2}} and \code{\link[geostan]{prep_car_data2}}. For more details, see: \code{vignette("raster-regression", package = "geostan")}.
#'
#' @return
#'
#' Returns the eigenvalues of the row-standardized spatial weights matrix (a numeric vector length `row * col`).
#'  
#' @seealso
#' \code{\link[geostan]{prep_sar_data2}}, \code{\link[geostan]{prep_car_data2}}
#' 
#' @source
#'
#' Griffith, Daniel A. (2000). Eigenfunction properties and approximations of selected incidence matrices employed in spatial analyses. *Linear Algebra and its Applications* 321 (1-3): 95-112. \doi{10.1016/S0024-3795(00)00031-8}.
#'
#' @examples
#'
#' e <- eigen_grid(row = 50, col = 95)
#' print(head(e, 25))
#'
#' @export
eigen_grid <- function(row = 5, col = 5) {
    P1 <- matrix(1, nrow = row)
    Q1 <- matrix(1, nrow = col)
    Pl <- eigen_1DW(P = row)
    Ql <- eigen_1DW(P = col)
    Ksum <- kronecker(P1, Ql) + kronecker(Pl, Q1)
    lambda <- sort( Ksum / 2 )
    return (lambda)
}

#' Prepare data for the CAR model: raster analysis
#'
#' @description Prepare a list of data required for the CAR model; this is for working with (large) raster data files only. For non-raster analysis, see \link[geostan]{prep_car_data}.
#'
#' @param row Number of rows in the raster 
#' @param col Number of columns in the raster
#' @param quiet Controls printing behavior. By default, `quiet = FALSE` and the range of permissible values for the spatial dependence parameter is printed to the console.
#' 
#' @details
#'
#' Prepare input data for the CAR model when your dataset consists of observations on a regular (rectangular) tessellation, such as a raster layer or remotely sensed imagery. The rook criteria is used to determine adjacency. This function uses Equation 5 from Griffith (2000) to generate approximate eigenvalues for a row-standardized spatial weights matrix from a P-by-Q dimension regular tessellation.
#'
#' This function can accommodate very large numbers of observations for use with \code{\link[geostan]{stan_car}}; for large N data, it is also recommended to use `slim = TRUE` or the `drop` argument. For more details, see: \code{vignette("raster-regression", package = "geostan")}.
#'
#'
#' @return A list containing all of the data elements required by the CAR model in \code{\link[geostan]{stan_car}}.
#' 
#' @source
#'
#' Griffith, Daniel A. (2000). Eigenfunction properties and approximations of selected incidence matrices employed in spatial analyses. *Linear Algebra and its Applications* 321 (1-3): 95-112. \doi{10.1016/S0024-3795(00)00031-8}. 
#'
#' @seealso
#' \code{\link[geostan]{prep_sar_data2}}, \code{\link[geostan]{prep_car_data}}, \code{\link[geostan]{stan_car}}.
#' 
#' @examples
#'
#' row = 100
#' col = 120
#' car_dl <- prep_car_data2(row = row, col = col)
#'
#' @export
#' @md
prep_car_data2 <- function(row = 100, col = 100, quiet = FALSE) {
    
    stopifnot(row > 2 & col > 2)
    N <- row * col
    Idx <- 1:N

    C_elements <- dims_to_W_elements(row = row, col = col)
    
    # create sparse matrix (row-standardized adjacency for WCAR)
    C <- Matrix::sparseMatrix(i = C_elements$I,
                              j = C_elements$J,
                              x = C_elements$x,
                              dims = c(N, N))

    # CRS representation
    car.dl <- list(A_w = C@x,
                A_v = C@i + 1,
                A_u = C@p + 1)

    # for wcar_normal_lpdf (including placeholders)
    car.dl$nA_w <- length(car.dl$A_w)
    car.dl$WCAR <- 1

    # Number of neighbors per row
    tmp <- aggregate(C_elements$x, by = list(I = C_elements$I), FUN = length)
    tmp <- tmp[sort(tmp$I), ]
    Ni <- tmp$x

    car.dl$Delta_inv <- Ni
    car.dl$log_det_Delta_inv <- sum(log(Ni))
    car.dl$n <- N

    # eigenvalues of the WCAR connectivity matrix
    car.dl$lambda <- eigen_grid(row = row, col = col)
    
    # limits of permissible/possible rho values
    car.dl$rho_lims <- 1 / range(car.dl$lambda)
    if (!quiet) {
            r_rho_lims <- round( car.dl$rho_lims, 3)
            message("Range of permissible rho values: ", r_rho_lims[1], ", ", r_rho_lims[2])
    }
    
    car.dl$style <- "WCAR"
    car.dl$C <- C
    return (car.dl)
}

#' Prepare data for SAR model: raster analysis
#'
#' @description Prepares a list of data required for using the SAR model; this is for working with (large) raster data files. For non-raster analysis, see \link[geostan]{prep_sar_data}.
#'
#' @param row Number of rows in the raster
#' @param col Number of columns in the raster
#' @param quiet Controls printing behavior. By default, `quiet = FALSE` and the range of permissible values for the spatial dependence parameter is printed to the console.
#'
#' @details
#'
#' Prepare data for the SAR model when your raw dataset consists of observations on a regular tessellation, such as a raster layer or remotely sensed imagery. The rook criteria is used to determine adjacency. This function uses Equation 5 from Griffith (2000) to calculate the eigenvalues for a row-standardized spatial weights matrix of a P-by-Q dimension regular tessellation.
#'
#' This function can accommodate very large numbers of observations for use with \code{\link[geostan]{stan_sar}}; for large N data, it is also recommended to use `slim = TRUE` or the `drop` argument. For details, see: \code{vignette("raster-regression", package = "geostan")}.
#'
#' @seealso
#' \code{\link[geostan]{prep_car_data2}}, \code{\link[geostan]{prep_sar_data}}, \code{\link[geostan]{stan_sar}}.
#'
#' @return A list containing all of the data elements required by the SAR model in \code{\link[geostan]{stan_sar}}.
#' 
#' @source
#'
#' Griffith, Daniel A. (2000). Eigenfunction properties and approximations of selected incidence matrices employed in spatial analyses. *Linear Algebra and its Applications* 321 (1-3): 95-112. \doi{10.1016/S0024-3795(00)00031-8}.
#'
#' @examples
#'
#' row = 100
#' col = 120
#' sar_dl <- prep_sar_data2(row = row, col = col)
#'
#' @export 
prep_sar_data2 <- function(row, col, quiet = FALSE) {
    stopifnot(row > 2 & col > 2)
    N <- row * col
    Idx <- 1:N

    W_elements <- dims_to_W_elements(row = row, col = col)
    
    # create sparse matrix: W
    W <- Matrix::sparseMatrix(i = W_elements$I,
                              j = W_elements$J,
                              x = W_elements$x,
                              dims = c(N, N))

    # CRS representation for Stan
    sar.dl <- rstan::extract_sparse_parts(W)
    names(sar.dl) <- paste0("W_", names(sar.dl))
    sar.dl$nW_w <- length(sar.dl$W_w)

    # eigenvalues of Wx
    sar.dl$eigenvalues_w <- eigen_grid(row = row, col = col)

    # N, rho limits
    sar.dl$n <- N    
    rho_lims <- 1/range(sar.dl$eigenvalues_w)
    if (!quiet) {
            r_rho_lims <- round( rho_lims, 3)
            message("Range of permissible rho values: ", r_rho_lims[1], ", ", r_rho_lims[2])
    }    
    sar.dl$rho_min <- min(rho_lims)
    sar.dl$rho_max <- max(rho_lims)
    sar.dl$W <- W
    return(sar.dl)
    }
   

#' Eigenvalues for W matrix for 1-dimensional connectivity, as in time series modeling.
#' 
#' The eigen_grid function (and thus prep_sar_data2, prep_car_data2) requires this function.
#' 
#' This implements Griffith 2000, Equation 4.
#'
#' @param P Number of observations.
#'
#' @details
#'
#' Observations are aligned in 1-dimensional space (a row, as in time series). Connect these by adjacency, construct the adjacency matrix, row-standardize the matrix, and then calculate eigenvalues of the matrix. This provides an approximation to the eigenvalues, and forms a part of the solution to the identification of the eigenvalues of two-dimensional regular tessellations.
#' 
#' @source
#'
#' Griffith, Daniel A. (2000). Eigenfunction properties and approximations of selected incidence matrices employed in spatial analyses. *Linear Algebra and its Applications* 321 (1-3): 95-112. \doi{10.1016/S0024-3795(00)00031-8}.
#' 
#' @noRd
eigen_1DW <- function(P = 5) {
        lambda <- NULL
    for (k in 1:( P - 1 )) {
        f = ( k * pi ) / ( P - 1 )
        cf <- cos( f ) 
        lambda <- c(lambda, cf)
    }
        lambda <- sort(c(lambda, 1))
        return (lambda)
}


#' given row and column dimensions, returns elements of row-standardized W matrix for regular (rectangular) tessellation with Rook adjacency
#'
#' This is required by prep_car_data2 and prep_sar_data2
#' @noRd
dims_to_W_elements <- function(row = 100, col = 100) {
    stopifnot(row > 2 & col > 2)
    N <- row * col
    Idx <- 1:N
    top_edge <- 2:(col - 1)
    bottom_edge <- (N-col+2):(N-1)
    right_edge <- col * 2:(row - 1)
    left_edge <- right_edge - (col - 1)
    top_left_corner <- 1
    top_right_corner <- col
    bottom_left_corner <- N - (col - 1)
    bottom_right_corner <- N
    
    # get row and column indices of all non-zero entries of W
    I_top_edge <- rep(top_edge, times = 3)
    J_top_edge <- c(top_edge - 1, top_edge + 1, top_edge + col)
    
    I_bottom_edge <-  rep(bottom_edge, time = 3)
    J_bottom_edge <- c(bottom_edge - 1, bottom_edge + 1, bottom_edge - col)
    
    I_right_edge <- rep(right_edge, times = 3)
    J_right_edge <- c(right_edge - 1, right_edge - col, right_edge + col)
    
    I_left_edge <- rep(left_edge, times = 3)
    J_left_edge <- c(left_edge + 1, left_edge - col, left_edge + col)
    
    I_top_left_corner <- rep(top_left_corner, times = 2)
    J_top_left_corner <- c(top_left_corner + 1, top_left_corner + col)
    
    I_bottom_left_corner <- rep(bottom_left_corner, times = 2)
    J_bottom_left_corner <- c(bottom_left_corner - col, bottom_left_corner + 1)
                                        
    I_top_right_corner <- rep(top_right_corner, times = 2)
    J_top_right_corner <- c(top_right_corner - 1, top_right_corner + col)

    I_bottom_right_corner <- rep(bottom_right_corner, times = 2)
    J_bottom_right_corner <- c(bottom_right_corner - 1, bottom_right_corner - col)

    I_edges <- c(
        I_top_edge,
        I_bottom_edge,
        I_right_edge,
        I_left_edge,
        I_top_left_corner,
        I_bottom_left_corner,
        I_top_right_corner,
        I_bottom_right_corner
    )

    J_edges <- c(
        J_top_edge,
        J_bottom_edge,
        J_right_edge,
        J_left_edge,
        J_top_left_corner,
        J_bottom_left_corner,
        J_top_right_corner,
        J_bottom_right_corner
    )

    middle <- Idx[-which(Idx %in% unique(I_edges))]
    I_middle <- rep(middle, times = 4)
    J_middle <- c(middle - 1, middle + 1, middle - col, middle + col)

    ##
    ## create W matrix
    I <- c(I_edges, I_middle)
    J <- c(J_edges, J_middle)
    
    # values for W
    x <- c(
        rep(1/3, times = length(c(I_top_edge,
                                I_bottom_edge,
                                I_left_edge,
                                I_right_edge))),
        rep(1/2, times = length(c(I_top_left_corner,
                                I_bottom_left_corner,
                                I_top_right_corner,
                                I_bottom_right_corner
                                ))),
        rep(1/4, times = length(I_middle))
    )

    # check sizes
    stopifnot(length(I) == length(J) & length(J) == length(x))

    W_elements <- data.frame(I = I,
                             J = J,
                             x = x)
    
    return (W_elements)
}

