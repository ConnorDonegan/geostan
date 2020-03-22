data {
  int<lower=0> n; // number of observations
  int<lower=0> dx; // number of covariates
  vector[n] y; // outcome variable
  matrix[n, dx] x; // covariates
  vector[3] alpha_prior; // other priors
  row_vector[dx] beta_prior[3];
  vector[3] sigma_prior;
  vector[3] alpha_tau_prior;
  vector[2] t_nu_prior;
  int<lower=0,upper=1> is_student; 
  int<lower=0,upper=1> has_re; // varying intercepts component
  int<lower=0> n_ids;
  int<lower=0,upper=n_ids> id[n];
}

transformed data {
  // use the QR decomposition on the matrix of covariates
  matrix[n, dx] Q_ast;
  matrix[dx, dx] R_ast;
  matrix[dx, dx] R_inverse;
  Q_ast = qr_Q(x)[, 1:dx] * sqrt(n - 1);
  R_ast = qr_R(x)[1:dx, ] / sqrt(n - 1);
  R_inverse = inverse(R_ast);
}

parameters {
  real<lower=0> sigma;
  real<lower=0> nu[is_student];
  real intercept;
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];
  vector[dx] beta_tilde;
}

transformed parameters {
  vector[dx] beta;
  vector[n] fitted;
  fitted = intercept + Q_ast * beta_tilde;
  beta = R_inverse * beta_tilde;
  if(has_re) {
    for (i in 1:n) {
      fitted[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
    }
  }
}

model {
  intercept ~ student_t(alpha_prior[1], alpha_prior[2], alpha_prior[3]);
  beta ~ student_t(beta_prior[1], beta_prior[2], beta_prior[3]);
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
  sigma ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]); 
  if (is_student) {
    nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);
    y ~ student_t(nu[1], fitted, sigma);
  } else {
    y ~ normal(fitted, sigma);
  }
}

generated quantities {
  vector[n] log_lik;
  vector[n] yrep;
  vector[n] residual;
  vector[n_ids] alpha_re;
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }
  for (i in 1:n) {
    residual[i] = y[i] - fitted[i];
    if (is_student) {
      log_lik[i] = student_t_lpdf(y[i] | nu[1], fitted[i], sigma);
      yrep[i] = student_t_rng(nu[1], fitted[i], sigma);
    } else {
      log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma);
      yrep[i] = normal_rng(fitted[i], sigma);
    }
  }
}

