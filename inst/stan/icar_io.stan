functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}

data {
  int<lower=0> n; // number of observations
  int<lower=0> n_edges;
  vector[n] y; // outcome variable
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  vector[3] alpha_prior; // other priors
  vector[3] sigma_prior;
  vector[3] alpha_tau_prior;
  vector[2] t_nu_prior;
  int<lower=0,upper=1> is_student; 
  int<lower=0,upper=1> has_re; // varying intercepts component
  int<lower=0> n_ids;
  int<lower=0,upper=n_ids> id[n];
}

parameters {
  real<lower=0> sigma;
  real<lower=0> nu[is_student];
  real intercept;
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];
  vector[n] phi;
}

transformed parameters {
  vector[n] fitted;
  fitted = intercept + phi;
  if(has_re) {
    for (i in 1:n) {
      fitted[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
    }
  }
}

model {
  phi ~ icar_normal_lpdf(n, node1, node2);
  intercept ~ student_t(alpha_prior[1], alpha_prior[2], alpha_prior[3]);
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

