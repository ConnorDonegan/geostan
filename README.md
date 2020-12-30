geostan
=======

Bayesian Spatial Analysis
-------------------------

The **geostan** R package provides a user-friendly interface to
hierarchical Bayesian models (HBMs) for areal data. It is designed for
ease of use with a particular emphasis on spatial epidemiology and survey data. All of the
models are fit using the Stan probabilistic programming language. The
package is still under development, you see more on the [package website](connordonegan.github.io/geostan)

The following models are implemented:

-   Generalized linear models with Gaussian, Student’s *t*, Poisson, and
    Binomial likelihood functions.
-   Conditional Autoregressive (CAR) models.
-   Eigenvector spatial filtering (ESF) models.
-   Intrinsic conditional autoregressive (IAR) models including Besag-York-Mollie (BYM) models and the innovation (BYM2) introduced by Riebler et al. (2016)

All of the models can incorporate additional ‘varying intercept’
terms for partial pooling of observations across geographies and users are encouraged to model data
uncertainty using, e.g., standard errors on estimates from the American Community Suvey.

### Package installation

If you’re using Linux you can try out a development version of the
package:

``` r
remotes::install_github("ConnorDonegan/geostan")
```


