# geostan 0.1.2

The measurement error models have been updated in three important respects:

  (1) There is now a prep_me_data function which must be used to create the list of data for the ME models. See `?prep_me_data`.
  (2) For covariates that are proportions or rates, the ME models now have an option for using a logit transformation on the variable. Again, see `?prep_me_data` for usage.
  (3) Previously, when using `stan_car`, ME models automatically employed the CAR model as a prior for the modeled covariates. That has changed, so that the default behavior for the ME models is the same across all `stan_*` models (CAR, GLM, ESF, ICAR). 

The second change addresses a limitation of the CAR prior models for the ME models. These are particularly important for variables that are highly skewed, such as the poverty rate. To determine whether a transformation should be considered, it can be helpful to evaluate results of the ME model (with the untransformed covariate) using the `me_diag` function. The logit transform is done on the 'latent' (modeled) variable, not the raw covariate. This transformation cannot be applied to the raw data by the user because the standard errors of covariate estimates (e.g., ACS standard errors) cannot be logit transformed.

# geostan 0.1.1

The stan files for the CAR model have been combined with the 'foundation.stan' file, which compresses the file size considerably. The vignette on spatial autocorrelation has also been updated to include model diagnostics, and a new example has been added to the stan_car documentation.

# geostan 0.1.0

geostan's first release.

