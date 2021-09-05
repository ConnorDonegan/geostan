// parameters for the process model //
  real intercept;
  vector[dwx] gamma;
  vector[dx_all] beta;
  real<lower=0> nu[is_student]; 
  real<lower=0> sigma[has_sigma];
// for partial pooling across groups/geographies
  vector[n_ids] alpha_re;
  real<lower=0> alpha_tau[has_re];
// observational error models //
  vector<lower=bounds[1],upper=bounds[2]>[n] x_true[dx_me];
  vector<lower=bounds[1],upper=bounds[2]>[dx_me] mu_x_true;
  vector<lower=0>[dx_me] sigma_x_true;
  vector<lower=ME_prior_car_rho[1], upper=ME_prior_car_rho[2]>[spatial_me ? dx_me : 0] car_rho_x_true;
  vector<lower=0>[spatial_me ? 0 : dx_me] nu_x_true;

