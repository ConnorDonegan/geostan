data {
#include parts/glm_data.stan
  int y[n]; // outcome variable
}

transformed data {
#include parts/QR.stan
}

parameters {
#include parts/glm_parameters.stan
}

transformed parameters {
#include parts/glm_trans_params.stan
}

model {
#include parts/glm_model.stan
  y ~ poisson_log(f); 
 }

generated quantities {
  vector[n] yrep;
  vector[n] fitted;
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
    log_lik[i] = poisson_log_lpmf(y[i] | f[i]);
  }
}

