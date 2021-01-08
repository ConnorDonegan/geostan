# geostan

The **geostan** R package provides a user-friendly interface to
hierarchical Bayesian models (HBMs) for spatial and spatio-temporal data. It is designed for
relative ease of use and developed with a particular emphasis on spatial epidemiology and survey data for public health.
All of the models are fit using the Stan probabilistic programming language, but users only need to be familiar with the R language.

The package and webiste are still under development but will be completed soon.

## Features

### Spatial data analysis tools

Tools for visualizing and measuring spatial autocorrelation and map patterns, for exploratory analysis and model diagnostics.

### Spatial regression 

Model count outcomes with a variety of conditional autoregressive models (CAR, ICAR, BYM, BYM2); model continuous or count outcomes with the eigenvector spatial filtering (ESF) methodology.

### Disease mapping and demography

Model small-area incidence rates or rate-ratios with mortality or disease data.

### Observational error models 

Incorporate critical information on data quality into any **geostan** model. Purpose-built for social, economic, and health surveys such as the American Community Survey (ACS).

### The RStan ecosystem

Full integration with a suite of state-of-the-art R packages for Bayesian inference and model exploration.
