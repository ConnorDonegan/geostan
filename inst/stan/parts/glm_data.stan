// design matrix, exchangeable random effects, and priors
  int<lower=0,upper=1> is_student;
  int<lower=0,upper=1> has_sigma;
  int<lower=0> n; // number of observations
  int<lower=0> dx; // number of covariates
  matrix[n, dx] x; // covariates
  vector[n] offset; // is a vector of zeros by default, otherwise a user-provided vector
  int<lower=0,upper=1> has_re; // has random effects? (or varying intercept)
  int<lower=0> n_ids; // number of random effects
  int<lower=0,upper=n_ids> id[n]; // identifier for the observational units associated with the random effects term
  vector[2] alpha_prior; // prior on the intercept
  row_vector[dx] beta_prior[2]; // coefficient priors
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts
  vector[2] t_nu_prior;
  vector[3] sigma_prior;

