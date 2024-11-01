  //ICAR
  vector[m] alpha_phi;  
  vector[type ? n : 0] phi_tilde;
  array[type ? 1 : 0] real<lower=0> spatial_scale;
  vector[type > 1 ? n : 0] theta_tilde;
  array[type == 2] real<lower=0> theta_scale;
  array[type == 3] real<lower=0,upper=1> rho;

  // ESF
  array[dev ? 1 : 0] real<lower=0> aux1_global;
  array[dev ? 1 : 0] real<lower=0> aux2_global;
  vector<lower=0>[dev] aux1_local;
  vector<lower=0>[dev] aux2_local;
  array[dev ? 1 : 0] real<lower=0> caux;
  vector[dev] z;

  // CAR/SAR
  vector[(car > 0 || sar > 0) && !is_auto_gaussian ? n : 0] log_lambda;
  array[car ? 1 : 0] real<lower=0> car_scale;
  array[car ? 1 : 0] real<lower=car_rho_lims[1], upper=car_rho_lims[2]> car_rho;     
  array[sar ? 1 : 0] real<lower=0> sar_scale;
  array[sar ? 1 : 0] real<lower=sar_rho_lims[1], upper=sar_rho_lims[2]> sar_rho;

  // GLM
  // parameters for the process model //
  real intercept;
  vector[dwx] gamma_qr;
  vector[dx_all] beta_qr;
  array[is_student] real<lower=0> nu; 
  array[has_sigma] real<lower=0> sigma;

  // for partial pooling across groups/geographies
  vector[n_ids] alpha_re;
  array[has_re] real<lower=0> alpha_tau;

  // observational error models (error in X) //
  array[dx_me] vector<lower=bounds[1],upper=bounds[2]>[n] x_true;
  vector[dx_me] mu_x_true;
  vector<lower=0>[dx_me] sigma_x_true;
  vector<lower=prior_rhox_true[1], upper=prior_rhox_true[2]>[spatial_me ? dx_me : 0] car_rho_x_true;
  vector<lower=0>[spatial_me ? 0 : dx_me] nu_x_true;

