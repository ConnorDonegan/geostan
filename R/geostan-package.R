#' The 'geostan' package.
#'
#' @description Bayesian spatial modeling powered by Stan. `geostan` offers access to a variety of hierarchical spatial models using the R formula interface. It is designed for spatial epidemiology and public health research but has broad applicability. 
#'
#' @docType package
#' @name geostan-package
#' @aliases geostan
#' @useDynLib geostan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling 
#'
#' @references
#'
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian Estimation of Spatial Filters with Moran’s Eigenvectors and Hierarchical Shrinkage Priors. Spatial Statistics. https://doi.org/10.1016/j.spasta.2020.100450
#'
#' Joseph, Max (2016). Exact Sparse CAR Models in Stan. Stan Case Studies, Vol. 3. https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. Spatial and spatio-temporal epidemiology, 31, 100301.
#' 
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL
