# geostan 0.3.0

## New additions

 - New exploratory spatial data analysis functions have been added: the Geary Ratio (GR) and the local Geary's C. These complement the Moran coefficient and local Moran's I. 
 - The vignette on spatial autocorrelation has been updated and expanded, including with a short discussion of exploratory spatial data analysis (ESDA). 
 - The vignette on spatial measurement error models/working with ACS data has been completely re-written.

## Minor changes

 - geostan models can now be used with the bridgesampling package for model camparison with Bayes factors (e.g., use `bridge_sampler(geostan_fit$stanfit)`). By default, geostan only collects MCMC samples for parameters that are expected to be of some interest for users. To become compatible with bridgesampling, the `keep_all` argument was added to all of the model fitting functions. For important background and details see the bridgesampling package documentation and vignettes on [CRAN](https://CRAN.R-project.org/package=bridgesampling).
 - stan_car now has an option to provide the connectivity matrix C, which is used to calculate spatial-lag of X (SLX) terms and residual spatial autocorrelation. Previously, there was no option to provide this matrix, as it was taken from the car_parts argument. However, that choice is only appropriate when the WCAR specification is used. Now, if C is missing and the WCAR specification has not been used a warning will appear.
 - Previously, the `lisa` function would automatically center and scale the variate before computing local Moran's I. Now, the variate will be centered and scaled by default but the user has the option to turn the scaling off (so the variate will be centered, but not divided by its standard deviation). This function also row-standardized the spatial weights matrix automatically, but there was no reason why. That's not done anymore.

# geostan 0.2.1

## Minor changes

The distance-based CAR models that are prepared by the `prep_car_data` function have changed slightly. The conditional variances were previously a function of the sum of neighboring inverse distances (in keeping with the specification of the connectivity matrix); this can lead to very skewed frequency distributions of the conditional variances. Now, the conditional variances are equal to the inverse of the number of neighboring sites. This is in keeping with the more common CAR model specifications.

# geostan 0.2.0

## Major changes

### Models for censored disease and mortality data

geostan now supports Poisson models with censored count data, a common problem in public health research where small area disease and mortality counts are censored below a threshold value. Model for censored outcome data can now be implemented using the `censor_point` argument found in all of the model fitting functions (stan_glm, stan_car, stan_esf, stan_icar).

### Measurement error models improved

The measurement error models have been updated in three important respects:

  - There is now a prep_me_data function which must be used to create the list of data for the ME models. See `?prep_me_data`.
  - For covariates that are proportions or rates, the ME models now have an option for using a logit transformation on the variable. Again, see `?prep_me_data` for usage.
  - Previously, when using `stan_car`, ME models automatically employed the CAR model as a prior for the modeled covariates. That has changed, so that the default behavior for the ME models is the same across all `stan_*` models (CAR, GLM, ESF, ICAR). 

The second change listed above is particularly useful for variables that are highly skewed, such as the poverty rate. To determine whether a transformation should be considered, it can be helpful to evaluate results of the ME model (with the untransformed covariate) using the `me_diag` function. The logit transform is done on the 'latent' (modeled) variable, not the raw covariate. This transformation cannot be applied to the raw data by the user because that would require the standard errors of covariate estimates (e.g., ACS standard errors) to be adjusted for the transformation.

## Minor changes

### A predict method for marginal effects

A `predict` method has been introduced for fitted geostan models; this is designed for calculating marginal effects. Fitted values of the model are still returned using `fitted` and the posterior predictive distribution is still accessible via `posterior_predict`.

### Centering covariates with measurement error models

The `centerx` argument has been updated to handle measurement error models for covariates. The centering now happens inside the Stan model so that the means of the modeled covariates (latent variables) are used instead of the raw data mean. 

# geostan 0.1.1

## Minor changes

  - The stan files for the CAR model have been combined with the 'foundation.stan' file, which compresses the file size considerably.
  - The vignette on spatial autocorrelation has also been updated to include model diagnostics.
  - A new example has been added to the stan_car documentation.

# geostan 0.1.0

geostan's first release.

