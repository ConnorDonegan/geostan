% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan_esf.R
\name{stan_esf}
\alias{stan_esf}
\title{Spatial filtering}
\source{
Chun, Y., D. A. Griffith, M. Lee and P. Sinha (2016). Eigenvector selection with stepwise regression techniques to construct eigenvector spatial filters. \emph{Journal of Geographical Systems}, 18(1), 67-85. \doi{10.1007/s10109-015-0225-3}.

Dray, S., P. Legendre & P. R. Peres-Neto (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). \emph{Ecological Modeling}, 196(3-4), 483-493.

Donegan, C., Y. Chun and A. E. Hughes (2020). Bayesian estimation of spatial filters with Moran’s Eigenvectors and hierarchical shrinkage priors. \emph{Spatial Statistics}. \doi{10.1016/j.spasta.2020.100450} (open access: \doi{10.31219/osf.io/fah3z}).

Griffith, Daniel A., and P. R. Peres-Neto (2006). Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses. \emph{Ecology} 87(10), 2603-2613.

Griffith, D., and Y. Chun (2014). Spatial autocorrelation and spatial filtering, Handbook of Regional Science. Fischer, MM and Nijkamp, P. eds.

Griffith, D., Chun, Y. and Li, B. (2019). \emph{Spatial Regression Analysis Using Eigenvector Spatial Filtering}. Elsevier.

Piironen, J and A. Vehtari (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. In \emph{Electronic Journal of Statistics}, 11(2):5018-5051.
}
\usage{
stan_esf(
  formula,
  slx,
  re,
  data,
  C,
  EV = make_EV(C, nsa = nsa, threshold = threshold),
  nsa = FALSE,
  threshold = 0.25,
  family = gaussian(),
  prior = NULL,
  ME = NULL,
  centerx = FALSE,
  censor_point,
  prior_only = FALSE,
  chains = 4,
  iter = 2000,
  refresh = 500,
  keep_all = FALSE,
  slim = FALSE,
  drop = NULL,
  pars = NULL,
  control = NULL,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{formula}{A model formula, following the R \code{\link[stats]{formula}} syntax. Binomial models are specified by setting the left hand side of the equation to a data frame of successes and failures, as in \code{cbind(successes, failures) ~ x}.}

\item{slx}{Formula to specify any spatially-lagged covariates. As in, \code{~ x1 + x2} (the intercept term will be removed internally). When setting priors for \code{beta}, remember to include priors for any SLX terms.}

\item{re}{To include a varying intercept (or "random effects") term, \code{alpha_re}, specify the grouping variable here using formula syntax, as in \code{~ ID}. Then, \code{alpha_re} is a vector of parameters added to the linear predictor of the model, and:

\if{html}{\out{<div class="sourceCode">}}\preformatted{alpha_re ~ N(0, alpha_tau)
alpha_tau ~ Student_t(d.f., location, scale).
}\if{html}{\out{</div>}}}

\item{data}{A \code{data.frame} or an object coercible to a data frame by \code{as.data.frame} containing the model data.}

\item{C}{Spatial connectivity matrix which will be used to calculate eigenvectors, if \code{EV} is not provided by the user. Typically, the binary connectivity matrix is best for calculating eigenvectors (i.e., using \code{C = shape2mat(shape, style = "B")}). This matrix will also be used to calculate residual spatial autocorrelation and any user specified \code{slx} terms; it will be row-standardized before calculating \code{slx} terms. See \code{\link[geostan]{shape2mat}}.}

\item{EV}{A matrix of eigenvectors from any (transformed) connectivity matrix, presumably spatial (see \code{\link[geostan]{make_EV}}). If \code{EV} is provided, still also provide a spatial weights matrix \code{C} for other purposes; \code{threshold} and \code{nsa} are ignored for user provided \code{EV}.}

\item{nsa}{Include eigenvectors representing negative spatial autocorrelation? Defaults to \code{nsa = FALSE}. This is ignored if \code{EV} is provided.}

\item{threshold}{Eigenvectors with standardized Moran coefficient values below this \code{threshold} value will be excluded from the candidate set of eigenvectors, \code{EV}. This defaults to \code{threshold = 0.25}, and is ignored if \code{EV} is provided.}

\item{family}{The likelihood function for the outcome variable. Current options are \code{family = gaussian()}, \code{student_t()} and \code{poisson(link = "log")}, and \code{binomial(link = "logit")}.}

\item{prior}{A named list of parameters for prior distributions (see \code{\link[geostan]{priors}}):  \describe{

\item{intercept}{The intercept is assigned a Gaussian prior distribution (see \code{\link[geostan]{normal}}}.

\item{beta}{Regression coefficients are assigned Gaussian prior distributions. Variables must follow their order of appearance in the model \code{formula}. Note that if you also use \code{slx} terms (spatially lagged covariates), and you use custom priors for \code{beta}, then you have to provide priors for the slx terms. Since slx terms are \emph{prepended} to the design matrix, the prior for the slx term will be listed first.
}

\item{sigma}{For \code{family = gaussian()} and \code{family = student_t()} models, the scale parameter, \code{sigma}, is assigned a (half-) Student's t prior distribution. The half-Student's t prior for \code{sigma} is constrained to be positive.}

\item{nu}{\code{nu} is the degrees of freedom parameter in the Student's t likelihood (only used when \code{family = student_t()}). \code{nu} is assigned a gamma prior distribution. The default prior is \code{prior = list(nu = gamma2(alpha = 3, beta = 0.2))}. }

\item{tau}{The scale parameter for random effects, or varying intercepts, terms. This scale parameter, \code{tau}, is assigned a half-Student's t prior. To set this, use, e.g., \code{prior = list(tau = student_t(df = 20, location = 0, scale = 20))}.}

\item{beta_ev}{The eigenvector coefficients are assigned the horseshoe prior (Piironen and Vehtari, 2017), parameterized by \code{global_scale} (to control overall prior sparsity), plus the degrees of freedom and scale of a Student's t model for any large coefficients (see \code{\link[geostan]{priors}}). To allow the spatial filter to account for a greater amount of spatial autocorrelation (i.e., if you find the residuals contain spatial autocorrelation), increase the global scale parameter (to a maximum of \code{global_scale = 1}).}
}}

\item{ME}{To model observational uncertainty (i.e. measurement or sampling error) in any or all of the covariates, provide a list of data as constructed by the \code{\link[geostan]{prep_me_data}} function.}

\item{centerx}{To center predictors on their mean values, use \code{centerx = TRUE}. If the ME argument is used, the modeled covariate (i.e., latent variable), rather than the raw observations, will be centered. When using the ME argument, this is the recommended method for centering the covariates.}

\item{censor_point}{Integer value indicating the maximum censored value; this argument is for modeling censored (suppressed) outcome data, typically disease case counts or deaths. For example, the US Centers for Disease Control and Prevention censors (does not report) death counts that are nine or fewer, so if you're using CDC WONDER mortality data you could provide \code{censor_point = 9}.}

\item{prior_only}{Draw samples from the prior distributions of parameters only.}

\item{chains}{Number of MCMC chains to estimate. Default \code{chains = 4}.}

\item{iter}{Number of samples per chain. Default \code{iter = 2000}.}

\item{refresh}{Stan will print the progress of the sampler every \code{refresh} number of samples. Defaults to \code{500}; set \code{refresh=0} to silence this.}

\item{keep_all}{If \code{keep_all = TRUE} then samples for all parameters in the Stan model will be kept; this is necessary if you want to do model comparison with Bayes factors and the \code{bridgesampling} package.}

\item{slim}{If \code{slim = TRUE}, then the Stan model will not collect the most memory-intensive parameters (including n-length vectors of fitted values, log-likelihoods, and ME-modeled covariate values). This will disable many convenience functions that are otherwise available for fitted \code{geostan} models, such as the extraction of residuals, fitted values, and spatial trends, WAIC, and spatial diagnostics, and ME diagnostics; many quantities of interest, such as fitted values and spatial trends, can still be calculated manually using given parameter estimates. The "slim" option is useful for data-intensive routines, such as regression with raster data, Monte Carlo studies, and measurement error models. For more control over which parameters are kept or dropped, use the \code{drop} argument instead of \code{slim}.}

\item{drop}{Provide a vector of character strings to specify the names of any parameters that you do not want MCMC samples for. Dropping parameters in this way can improve sampling speed and reduce memory usage. The following parameter vectors can potentially be dropped from ESF models:
\describe{
\item{fitted}{The N-length vector of fitted values}
\item{alpha_re}{Vector of 'random effects'/varying intercepts.}
\item{x_true}{N-length vector of 'latent'/modeled covariate values created for measurement error (ME) models.}
\item{esf}{The N-length eigenvector spatial filter.}
\item{beta_ev}{The vector of coefficients for the eigenvectors.}
}
If \code{slim = TRUE}, then \code{drop} will be ignored.}

\item{pars}{Optional; specify any additional parameters you'd like stored from the Stan model.}

\item{control}{A named list of parameters to control the sampler's behavior. See \link[rstan]{stan} for details.}

\item{quiet}{By default, any prior distributions that have not been assigned by the user are printed to the console. If \code{quiet = TRUE}, these will not be printed.}

\item{...}{Other arguments passed to \link[rstan]{sampling}.}
}
\value{
An object of class class \code{geostan_fit} (a list) containing:
\describe{
\item{summary}{Summaries of the main parameters of interest; a data frame}
\item{diagnostic}{Widely Applicable Information Criteria (WAIC) with a measure of effective number of parameters (\code{eff_pars}) and mean log pointwise predictive density (\code{lpd}), and mean residual spatial autocorrelation as measured by the Moran coefficient.}
\item{data}{a data frame containing the model data}
\item{EV}{A matrix of eigenvectors created with \code{w} and \code{geostan::make_EV}}
\item{C}{The spatial weights matrix used to construct EV}
\item{family}{the user-provided or default \code{family} argument used to fit the model}
\item{formula}{The model formula provided by the user (not including ESF component)}
\item{slx}{The \code{slx} formula}
\item{re}{A list containing \code{re},  the random effects (varying intercepts) formula if provided, and
\code{data} a data frame with columns \code{id}, the grouping variable, and \code{idx}, the index values assigned to each group.}
\item{priors}{Prior specifications.}
\item{x_center}{If covariates are centered internally (\code{centerx = TRUE}), then \code{x_center} is a numeric vector of the values on which covariates were centered.}
\item{ME}{The \code{ME} data list, if one was provided by the user for measurement error models.}
\item{spatial}{A data frame with the name of the spatial component parameter ("esf") and method ("ESF")}
\item{stanfit}{an object of class \code{stanfit} returned by \code{rstan::stan}}
}
}
\description{
Fit a spatial regression model using eigenvector spatial filtering (ESF).
}
\details{
Eigenvector spatial filtering (ESF) is a method for spatial regression analysis. ESF is extensively covered in Griffith et al. (2019). This function implements the methodology introduced in Donegan et al. (2020), which uses Piironen and Vehtari's (2017) regularized horseshoe prior.

ESF decomposes spatial autocorrelation into a linear combination of various patterns, typically at different scales (such as local, regional, and global trends). By adding a spatial filter to a regression model, these spatial autocorrelation patterns are shifted from the residuals to the spatial filter. ESF models take the spectral decomposition of a transformed spatial connectivity matrix, \eqn{C}. The resulting eigenvectors, \eqn{E}, are mutually orthogonal and uncorrelated map patterns. The spatial filter equals \eqn{E \beta_{E}} where \eqn{\beta_{E}} is a vector of coefficients.

ESF decomposes the data into a global mean, \eqn{\alpha}, global patterns contributed by covariates \eqn{X \beta}, spatial trends \eqn{E \beta_{E}}, and residual variation. Thus, for \code{family=gaussian()},

\deqn{
y \sim Gauss(\alpha + X * \beta + E \beta_{E}, \sigma).
}

An ESF component can be incorporated into the linear predictor of any generalized linear model. For example, a spatial Poisson model for rare disease incidence may be specified as follows:

\deqn{y \sim Poisson(e^{O + \mu})}
\deqn{\mu = \alpha + E \beta_{E} + A }
\deqn{ A \sim Guass(0, \tau) }
\deqn{ \tau \sim student(20, 0, 2) }
\deqn{ \beta_{E} \sim horseshoe(.) }

The form of this model is similar to the BYM model (see \link[geostan]{stan_icar}), in the sense that it contains a spatially structured trend term (\eqn{E \beta_{E}}) and an unstructured ('random effects') term (\eqn{A}).

The \link[geostan]{spatial.geostan_fit} method will return \eqn{E \beta_{E}}.

The model can also be extended to the space-time domain; see \link[geostan]{shape2mat} to specify a space-time connectivity matrix.

The coefficients \eqn{\beta_{E}} are assigned the regularized horseshoe prior (Piironen and Vehtari, 2017), resulting in a relatively sparse model specification. In addition, numerous eigenvectors are automatically dropped because they represent trace amounts of spatial autocorrelation (this is controlled by the \code{threshold} argument). By default, \code{stan_esf} will drop all eigenvectors representing negative spatial autocorrelation patterns. You can change this behavior using the \code{nsa} argument.
\subsection{Additional functionality}{

The CAR models can also incorporate spatially-lagged covariates, measurement/sampling error in covariates (particularly when using small area survey estimates as covariates), missing outcome data, and censored outcomes (such as arise when a disease surveillance system suppresses data for privacy reasons). For details on these options, please see the Details section in the documentation for \link[geostan]{stan_glm}.
}
}
\examples{
\donttest{
data(sentencing)
# spatial weights matrix with binary coding scheme
C <- shape2mat(sentencing, style = "B", quiet = TRUE)

# expected number of sentences
log_e <- log(sentencing$expected_sents)

# fit spatial Poisson model with ESF + unstructured 'random effects'
fit.esf <- stan_esf(sents ~ offset(log_e),
                   re = ~ name,
                   family = poisson(),
                   data = sentencing,
                   C = C, 
                   chains = 2, iter = 800) # for speed only

# spatial diagnostics 
sp_diag(fit.esf, sentencing)

# plot marginal posterior distributions of beta_ev (eigenvector coefficients)
plot(fit.esf, pars = "beta_ev")

# calculate log-standardized incidence ratios (SIR)
#  # SIR = observed/exected number of cases
# in this case, prison sentences
library(ggplot2)
library(sf)
f <- fitted(fit.esf, rates = FALSE)$mean
SSR <-  f / sentencing$expected_sents
log.SSR <- log( SSR, base = 2 )

# map the log-SSRs
 ggplot(sentencing) +
   geom_sf(aes(fill = log.SSR)) +
   scale_fill_gradient2(
    midpoint = 0,
    name = NULL,
    breaks = seq(-3, 3, by = 0.5)
  ) +
   labs(title = "Log-Standardized Sentencing Ratios",
    subtitle = "log( Fitted/Expected ), base 2"
 ) +
   theme_void()
}
}
\author{
Connor Donegan, \email{connor.donegan@gmail.com}
}
