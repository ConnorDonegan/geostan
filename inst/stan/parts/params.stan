// parameters for the process model //
// regression model
  real intercept;
  vector[dwx] gamma;
  vector[dx_all] beta;
  real<lower=0> nu[is_student]; 
  real<lower=0> sigma[has_sigma];
// for partial pooling across groups/geographies
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];
// observational error models //
// covariates to model: bounded
  matrix<lower=bounds[1],upper=bounds[2]>[n, dx_me_bounded] x_true_bounded;
  vector<lower=bounds[1],upper=bounds[2]>[dx_me_bounded ? dx_me_bounded : 0] mu_x_true_bounded;
  vector<lower=0>[dx_me_bounded ? dx_me_bounded : 0] sigma_x_true_bounded;
  vector<lower=0>[dx_me_bounded ? dx_me_bounded : 0] nu_x_true_bounded;  
// covariates to model: unbounded
  matrix[n, dx_me_unbounded] x_true_unbounded;
  vector[dx_me_unbounded ? dx_me_unbounded : 0] mu_x_true_unbounded;
  vector<lower=0>[dx_me_unbounded ? dx_me_unbounded : 0] sigma_x_true_unbounded;
  vector<lower=0>[dx_me_unbounded ? dx_me_unbounded : 0] nu_x_true_unbounded;  
// offset model
  vector<lower=0>[model_offset ? n : 0] offset_est;
  real<lower=0> mu_offset[model_offset ? 1: 0];  
  real<lower=0> sigma_offset[model_offset ? 1 : 0];
  real<lower=0> nu_offset[model_offset ? 1 : 0];  
