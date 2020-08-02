// covariates to model, bounded and unbounded 
  matrix<lower=bounds[1],upper=bounds[2]>[n, dx_me_bounded] x_true_bounded;
  matrix[n, dx_me_unbounded] x_true_unbounded;
  vector<lower=0>[model_offset ? n : 0] offset_est;  
// regression parameters
  real intercept;
  vector[dwx] gamma;
  vector[dx_all] beta;
  real<lower=0> nu[is_student]; 
  real<lower=0> sigma[has_sigma];
// for partial pooling across groups/geographies
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];


