# geostan <img src="man/figures/logo.png" align="right" width="160" />

The **geostan** R package supports a complete spatial analysis
workflow with hierarchical Bayesian models (HBMs) for areal
data and a variety of functions for visualizing spatial data and model results.

### Disease mapping and spatial regression

Model small-area incidence rates with mortality or disease data recorded across areal units like counties or census tracts.

### Observational uncertainty 

Incorporate information on data reliability into any **geostan** model. Built specifically for American Community Survey (ACS) data.

### Spatial analysis tools

Tools for visualizing and measuring spatial autocorrelation and map patterns, for exploratory analysis and model criticism.

### Custom Stan models

Tools for building custom spatial models in [Stan](https://mc-stan.org/).

### RStan ecosystem

Compatible with a suite of high-quality R packages for Bayesian inference.

## Installation

Install **geostan** from CRAN using:

``` r
install.packages("geostan")
```

## Citation

Please cite **geostan** if you use it in your work, and cite any appropriate methodology papers as well (found in the package documentation page for the functions you use).

 * Donegan, Connor (2021). geostan: Bayesian Spatial Analysis. R package Version 0.1.1 https://connordonegan.github.io/geostan
 
 * Donegan, Connor, Yongwan Chun, and Daniel A. Griffith. Modeling community health with areal data: Bayesian inference with survey standard errors and spatial structure. International Journal of Environmental Research and Public Health 18.13 (2021): 6856. DOI: 10.3390/ijerph18136856

All **geostan** models are built using **Stan**, so be sure to cite **Stan** too:

 * Carpenter B., Gelman A., Hoffman M. D., Lee D., Goodrich B., Betancourt M., Brubaker M., Guo J., Li P., and Riddell A. (2017). Stan: A probabilistic programming language. Journal of Statistical Software. 76(1). DOI: 10.18637/jss.v076.i01

<br />

<span style="color:gray">The **geostan** package was built with the help of **rstantools**:</span>

<span style="color:gray">  Gabry, Jonah, Ben Goodrich, and Martin Lysy (2021). rstantools: Tools for Developing R Packages Interfacing with 'Stan'. R package version 2.1.1 https://mc-stan.org/rstantools/index.html</span>.
