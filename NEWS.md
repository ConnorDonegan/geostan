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

