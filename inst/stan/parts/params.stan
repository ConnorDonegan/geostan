  //ICAR
  vector[m] alpha_phi;  
  vector[type ? n : 0] phi_tilde;
  real<lower=0> spatial_scale[type ? 1 : 0];
  vector[type > 1 ? n : 0] theta_tilde;
  real<lower=0> theta_scale[type == 2];
  real<lower=0,upper=1> rho[type == 3];
  // ESF
  real<lower=0> aux1_global[dev ? 1 : 0];
  real<lower=0> aux2_global[dev ? 1 : 0];
  vector<lower=0>[dev] aux1_local;
  vector<lower=0>[dev] aux2_local;
  real<lower=0> caux[dev ? 1 : 0];
  vector[dev] z;
  // CAR
  vector[car && !is_auto_gaussian ? n : 0] log_lambda;
  real<lower=0> car_scale[car ? 1 : 0];
  real<lower=car_rho_lims[1], upper=car_rho_lims[2]> car_rho[car ? 1 : 0];     
  // GLM
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
  vector<lower=prior_rhox_true[1], upper=prior_rhox_true[2]>[spatial_me ? dx_me : 0] car_rho_x_true;
  vector<lower=0>[spatial_me ? 0 : dx_me] nu_x_true;

