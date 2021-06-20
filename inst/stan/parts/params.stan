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
  vector<lower=bounds[1],upper=bounds[2]>[n] x_true_bounded[dx_me_bounded];
  vector<lower=bounds[1],upper=bounds[2]>[dx_me_bounded] mu_x_true_bounded;
  vector<lower=0>[dx_me_bounded] sigma_x_true_bounded;
  vector<lower=1/min(lambda), upper=1/max(lambda)>[spatial_me ? dx_me_bounded : 0] car_rho_x_true_bounded;
  vector<lower=0>[spatial_me ? 0 : dx_me_bounded] nu_x_true_bounded;
// covariates to model: unbounded
  vector[n] x_true_unbounded[dx_me_unbounded]; 
  vector[dx_me_unbounded] mu_x_true_unbounded;
  vector<lower=0>[dx_me_unbounded] sigma_x_true_unbounded;
  vector<lower=1/min(lambda), upper=1/max(lambda)>[spatial_me ? dx_me_unbounded : 0] car_rho_x_true_unbounded;
  vector<lower=0>[spatial_me ? 0 : dx_me_unbounded] nu_x_true_unbounded;

