#' The geostan R package.
#'
#' @description Bayesian spatial modeling powered by Stan. \code{geostan} offers access to a variety of hierarchical spatial models using the R formula interface. It is designed primarily for spatial epidemiology and public health research but is generally applicable to modeling areal data. Unique features of the package include its suitability for inference with small area estimates from the American Community Survey (using spatial observational error models), and its eigenvector spatial filtering methodology.
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
#' Carpenter, B., Gelman, A., Hoffman, M.D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., Riddell, A., 2017. Stan: A probabilistic programming language. Journal of statistical software 76. \doi{10.18637/jss.v076.i01}.
#'
#' Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. *Spatial Statistics*. \doi{10.1016/j.spasta.2020.100450}.
#'
#' Donegan, Connor and Chun, Yongwan and Griffith, Daniel A. (2021). Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. *Int. J. Env. Res. and Public Health* 18 (13): 6856. \doi{10.3390/ijerph18136856}. Data and code: \url{https://github.com/ConnorDonegan/survey-HBM}.
#' 
#' Gabry, J., Goodrich, B. and Lysy, M. (2020). rstantools: Tools for developers of R packages interfacing with Stan. R package version 2.1.1 \url{https://mc-stan.org/rstantools/}.
#' 
#' Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A., & DiMaggio, C. (2019). Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. Spatial and spatio-temporal epidemiology, 31, 100301. \doi{10.1016/j.sste.2019.100301}.
#' 
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. \url{https://mc-stan.org}
#'
NULL
