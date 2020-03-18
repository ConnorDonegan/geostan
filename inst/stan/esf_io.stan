data {
  int<lower=0> n; // number of observations
  int<lower=0> dev; // number of EVs
  vector[n] y; // outcome variable
  matrix[n, dev] EV; // eigenvectors
  real<lower=0> scale_global;  // horseshoe hyperparameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df; 
  vector[3] alpha_prior; // other priors
  vector[3] sigma_prior;
  vector[2] t_nu_prior;
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts
  int<lower=0,upper=1> is_student;
  int<lower=0,upper=1> has_re; // does the model have varying intercept?
  int<lower=0> n_ids; // number of 'regions' with varying intercpts
  int<lower=0,upper=n_ids> id[n]; // index for the regions
}

transformed data {
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  nu_global = 1;
  nu_local = 1;
}

parameters {
  real<lower=0> sigma;
  real<lower=0> nu[is_student];
  real intercept;
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[dev] aux1_local;
  vector<lower=0>[dev] aux2_local;
  real<lower=0> caux;
  vector[dev] z;
}

transformed parameters {
  // RHS prior on the EV matrix 
  real <lower=0> tau;
  vector<lower=0>[dev] lambda;
  vector<lower=0>[dev] lambda_tilde;
  vector[dev] beta_ev;
  real <lower=0> c;
  vector[n] f;
  tau = aux1_global * sqrt(aux2_global) * scale_global * sigma;
  lambda = aux1_local .* sqrt(aux2_local);
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
  beta_ev = z .* lambda_tilde * tau;
  f = intercept + EV * beta_ev;
  if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
    }
  }
}

model {
  z ~ normal(0, 1);
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  intercept ~ student_t(alpha_prior[1], alpha_prior[2], alpha_prior[3]);
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
  sigma ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]); 
  if (is_student) {
    nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);
    y ~ student_t(nu[1], f, sigma);
  } else {
    y ~ normal(f, sigma);
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] yrep;
  vector[n] residual;
  vector[n] fitted;
  vector[n] esf;
  vector[n_ids] alpha_re;
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }
  for (i in 1:n) {
    fitted[i] = f[i];
    residual[i] = y[i] - fitted[i];
    esf[i] = EV[i] * beta_ev;
    if (is_student) {
      log_lik[i] = student_t_lpdf(y[i] | nu[1], fitted[i], sigma);
      yrep[i] = student_t_rng(nu[1], fitted[i], sigma);
    } else {
      log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma);
      yrep[i] = normal_rng(fitted[i], sigma);
    }
  }
}

