data { 
#include parts/glm_data.stan
#include parts/esf_data.stan
#include parts/binomial_data.stan
}

transformed data {
#include parts/QR.stan
}

parameters {
#include parts/esf_parameters.stan
#include parts/glm_parameters.stan
}

transformed parameters {
  vector<lower=0, upper=1>[n] p;
#include parts/glm_trans_params.stan
  p = inv_logit(f);  
}


model {
#include parts/glm_model.stan
  y ~ binomial(N, p); 
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
    fitted[i] = p[i] * N[i];
    residual[i] = fitted[i] - y[i];
    yrep[i] = binomial_rng(N[i], p[i]);
    log_lik[i] = binomial_lpmf(y[i] | N[i], p[i]);
  }
}
