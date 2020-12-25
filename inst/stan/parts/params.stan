// parameters for the process model //
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
  vector<lower=bounds[1],upper=bounds[2]>[dx_me_bounded] mu_x_true_bounded;
  vector<lower=0>[dx_me_bounded] sigma_x_true_bounded;
  vector<lower=0>[dx_me_bounded] nu_x_true_bounded;    
// covariates to model: unbounded
  vector<lower=bounds[1],upper=bounds[2]>[dx_me_unbounded] mu_x_true_unbounded;
  vector<lower=0>[dx_me_unbounded] sigma_x_true_unbounded;
  vector<lower=0>[dx_me_unbounded] nu_x_true_unbounded;
  
// spatial covariate models
  vector[n] x_true_unbounded[dx_me_unbounded]; 
  vector<lower=bounds[1],upper=bounds[2]>[n] x_true_bounded[dx_me_bounded]; 
  // unbounded
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] aux1_global_me_unbounded;
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] aux2_global_me_unbounded;
  vector<lower=0>[dev] aux1_local_me_unbounded[spatial_me ? dx_me_unbounded : 0];
  vector<lower=0>[dev] aux2_local_me_unbounded[spatial_me ? dx_me_unbounded : 0];
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] caux_me_unbounded;
  vector[dev] z_ev_me_unbounded[spatial_me ? dx_me_unbounded : 0];
  // bounded
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] aux1_global_me_bounded;
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] aux2_global_me_bounded;
  vector<lower=0>[dev] aux1_local_me_bounded[spatial_me ? dx_me_bounded : 0];
  vector<lower=0>[dev] aux2_local_me_bounded[spatial_me ? dx_me_bounded : 0];
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] caux_me_bounded;
  vector[dev] z_ev_me_bounded[spatial_me ? dx_me_bounded : 0];
// offset model
  vector<lower=0>[model_offset ? n : 0] offset_est;
  real<lower=0> mu_offset[model_offset ? 1: 0];  
  real<lower=0> sigma_offset[model_offset ? 1 : 0];
  real<lower=0> nu_offset[model_offset ? 1 : 0];  

