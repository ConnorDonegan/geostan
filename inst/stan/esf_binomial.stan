functions {
#include parts/priors.stan
}

data { // binomial esf
#include parts/all_data.stan
#include parts/esf_data.stan
  int y[n]; // successes
  int N[n]; // trials
}

transformed data {
}

parameters {
#include parts/esf_parameters.stan
#include parts/glm_parameters.stan
}

transformed parameters {
  vector[n] f;
  vector<lower=0, upper=1>[n] p;
  vector[dev] beta_ev;
  beta_ev = rhs_prior(dev, z, aux1_global, aux2_global, aux1_local, aux2_local, caux, scale_global, slab_scale, 1);
  f = offset + intercept + EV * beta_ev;
  if (dx) f += x * beta;
  if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
   }
  }
  p = inv_logit(f);
}

model {
#include parts/glm_model.stan
#include parts/rhs_model.stan
  y ~ binomial(N, p); 
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
    fitted[i] = p[i] * N[i];
    residual[i] = fitted[i] - y[i];
    yrep[i] = binomial_rng(N[i], p[i]);
    esf[i] = EV[i]*beta_ev;
    log_lik[i] = binomial_lpmf(y[i] | N[i], p[i]);
  }
}

