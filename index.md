### geostan

The **geostan** R package provides a user-friendly interface to
hierarchical Bayesian models (HBMs) for spatial and spatio-temporal data. It is designed for
relative ease of use and developed with a particular emphasis on spatial epidemiology and survey data for public health.
All of the models are fit using the Stan probabilistic programming language, but users only need to be familiar with the R language.

The following models are implemented:

-   Generalized linear models with Gaussian, Student’s *t*, Poisson, and
    Binomial likelihood functions.
-   Conditional Autoregressive (CAR) models.
-   Eigenvector spatial filtering (ESF) models.
-   Intrinsic conditional autoregressive (IAR) models including Besag-York-Mollie (BYM) models and the innovation (BYM2) introduced by Riebler et al. (2016)

All of the models can incorporate additional ‘varying intercept’
terms for partial pooling of observations across geographies and users are encouraged to model observational error using, e.g., standard errors on estimates from the American Community Suvey.

### Package installation

This package is still under development. If you’re using Linux you can try out a development version of the
package:

``` r
remotes::install_github("ConnorDonegan/geostan")
```

