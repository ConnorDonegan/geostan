data {
  int<lower=0> n; // number of observations
  int<lower=0> dx; // number of covariates
  int<lower=0> dev; // number of eigenvectors 
  matrix[n, dev] EV; // eigenvectors
  matrix[n, dx] x; // covariates
  int y[n]; // outcome variable
  vector[n] offset; // is a vector of zeros by default, otherwise a user-provided vector
  int<lower=0,upper=1> has_re; // has random effects? (or varying intercept)
  int<lower=0> n_ids; // number of random effects
  int<lower=0,upper=n_ids> id[n]; // identifier for the observational units associated with the random effects term
  real<lower=0> scale_global;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df; 
  vector[3] alpha_prior; // prior on the intercept
  row_vector[dx] beta_prior[3];
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts
}

transformed data {
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  vector[n] log_E;
  nu_global = 1;
  nu_local = 1;
  log_E = log(offset);
}

parameters {
  real intercept;
  vector[dx] beta;
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
  tau = aux1_global * sqrt(aux2_global) * scale_global;
  lambda = aux1_local .* sqrt(aux2_local);
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
  beta_ev = z .* lambda_tilde * tau;
  f = log_E + intercept + EV * beta_ev;
  if (dx) f += x * beta;
  if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
   }
  }
}

model {
  y ~ poisson_log(f); 
  z ~ std_normal();
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  intercept ~ student_t(alpha_prior[1], alpha_prior[2], alpha_prior[3]);
  if (dx) beta ~ student_t(beta_prior[1], beta_prior[2], beta_prior[3]);
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
}

generated quantities {
  vector[n] yrep;
  vector[n] fitted;
  vector[n] esf;
  vector[n] residual;
  vector[n] log_lik;
  vector[n_ids] alpha_re;
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }
  for (i in 1:n) {
    fitted[i] = exp(f[i]);
    residual[i] = fitted[i] - y[i];
    yrep[i] = poisson_log_rng(f[i]);
    esf[i] = EV[i]*beta_ev;
    log_lik[i] = poisson_log_lpmf(y[i] | f[i]);
  }
}

