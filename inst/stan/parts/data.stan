// likelihood function
  int<lower=0,upper=4> family;
// number of observations
  int<lower=0> n; 
// lower, upper bound for bounded observational error models
  vector[2] bounds;
// offest with measurement error information   
  vector<lower=0>[n] offset_obs; 
  vector<lower=0>[n] offset_me;
  int<lower=0,upper=1> model_offset;
// connectivity matrix
  int<lower=0> dwx;
  int wx_idx[dwx ? dwx : 1];
  int<lower=0> dw_nonzero;
  matrix[dwx ? n : 1, dwx ? n : 1] W;
// covariates with measurement error information
  int<lower=0> dx_obs;
  int<lower=0> dx_me_bounded;
  int<lower=0> dx_me_unbounded;
  int<lower=0> x_obs_idx[dx_obs ? dx_obs : 1];
  int<lower=0> x_me_bounded_idx[dx_me_bounded ? dx_me_bounded : 1];
  int<lower=0> x_me_unbounded_idx[dx_me_unbounded ? dx_me_unbounded : 1];
  matrix[n, dx_obs ? dx_obs : 0] x_obs;
  matrix<lower=bounds[1],upper=bounds[2]>[n, dx_me_bounded ? dx_me_bounded : 1] x_me_bounded;
  matrix[n, dx_me_unbounded ? dx_me_unbounded : 1] x_me_unbounded;
  matrix<lower=0>[n, dx_me_bounded ? dx_me_bounded : 1] sigma_me_bounded;
  matrix<lower=0>[n, dx_me_unbounded ? dx_me_unbounded : 1] sigma_me_unbounded;
// exchangeable random effects
  int<lower=0,upper=1> has_re; // has random effects? (or varying intercept)
  int<lower=0> n_ids; // number of random effects
  int<lower=0,upper=n_ids> id[n]; // identifier for the observational units associated with the random effects term
// priors
  vector[2] alpha_prior; // prior on the intercept
  int<lower=0> dbeta_prior;  
  row_vector[dbeta_prior] beta_prior[2]; // coefficient priors, with any SLX terms listed first; 
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts
  vector[2] t_nu_prior;
  vector[3] sigma_prior;
  vector[n] y;
  int<lower=0> y_int[n];
  int<lower=0> trials[n];


