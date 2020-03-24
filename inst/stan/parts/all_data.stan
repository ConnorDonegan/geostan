// design matrix, exchangeable random effects, and priors
  int<lower=1,upper=4> family; // gauss=1,t=2,pois=3,binom=4
  int<lower=0> n; // number of observations
  int<lower=0> dx; // number of covariates
  matrix[n, dx] x; // covariates
  vector[n] offset; // is a vector of zeros by default, otherwise a user-provided vector
  int<lower=0,upper=1> has_re; // has random effects? (or varying intercept)
  int<lower=0> n_ids; // number of random effects
  int<lower=0,upper=n_ids> id[n]; // identifier for the observational units associated with the random effects term
  vector[3] alpha_prior; // prior on the intercept
  row_vector[dx] beta_prior[3]; // coefficient priors
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts

