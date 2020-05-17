functions {
#include parts/priors.stan
}

data {
#include parts/glm_data.stan
#include parts/icar_data.stan
  real<lower=0> scaling_factor; // scales the spatial component
  int y[n]; // outcome variable
}

transformed data {
#include parts/QR.stan
}

parameters {
#include parts/glm_parameters.stan
  vector[n] v;
  vector[n] u;
  real<lower=0> sigma_re;
  real logit_rho;
}

transformed parameters {
#include parts/bym2_trans_params.stan
}

model {
#include parts/glm_model.stan
#include parts/bym2_model.stan
  y ~ poisson_log(f); 
}

generated quantities {
  vector[n] yrep;
  vector[n] fitted;
  vector[n] residual;
  vector[n] log_lik;
  vector[n_ids] alpha_re;
  vector[n] phi; //scaled v
  vector[n] theta; // scaled u
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }
  for (i in 1:n) {
    phi[i] = sigma_re * sqrt(rho / scaling_factor) * v[i];
    theta[i] = sigma_re * sqrt(1 - rho) * u[i];
    fitted[i] = exp(f[i]);
    residual[i] = fitted[i] - y[i];
    log_lik[i] = poisson_log_lpmf(y[i] | f[i]);
    if (f[i] > 20) {
       print("f[i] too large (>20) for poisson_log_rng");
       yrep[i] = -1;
       } else {
    yrep[i] = poisson_log_rng(f[i]);
    }
  }
}



