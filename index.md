# geostan <img src="man/figures/logo.png" align="right" width="120" />

The **geostan** R package provides a user-friendly interface to
hierarchical Bayesian models (HBMs) for spatial and spatio-temporal data. It is designed for
relative ease of use and developed with a particular emphasis on spatial epidemiology and survey data for public health.
All of the models are built with RStan, but users only need to be familiar with the R language.

The package and webiste are still under development but will be completed soon.

## Features

### Spatial data analysis tools

Tools for visualizing and measuring spatial autocorrelation and map patterns, for exploratory analysis and model diagnostics.

### Spatial regression 

Model count outcomes with a variety of (intrinsic) conditional autoregressive models (CAR, ICAR, BYM, BYM2); model continuous or count outcomes with the eigenvector spatial filtering (ESF) methodology.

### Disease mapping and spatial demography

Model small-area incidence rates or standardized rate-ratios with mortality or disease data.

### Spatio-temporal models

Model spatial processes unfolding over time.

### Observational error models 

Incorporate information on data quality into any **geostan** model. Purpose-built for social, economic, and health surveys such as the American Community Survey (ACS).

### The RStan ecosystem

Compatible with a suite of high-quality R packages for Bayesian inference and model exploration.
