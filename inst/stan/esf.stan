functions {
#include parts/priors.stan
}
data {
#include parts/data.stan
// esf data
  int<lower=0> dev; // number of eigenvectors : now included in parts/data.stan
  matrix[n, dev] EV; // the eigenvectors : now included in parts/data.stan
  real<lower=0> scale_global;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
}

transformed data {
#include parts/trans_data.stan
}

parameters {
// ESF parameters
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[dev] aux1_local;
  vector<lower=0>[dev] aux2_local;
  real<lower=0> caux;
  vector[dev] z;
#include parts/params.stan
}

transformed parameters {
  vector[dev] beta_ev;  
  real error_scale;
#include parts/trans_params_declaration.stan
  if (has_sigma) {
    error_scale = sigma[1];
      } else {
    error_scale = 1;
      }
  beta_ev = rhs_prior(dev, z, aux1_global, aux2_global, aux1_local, aux2_local, caux, scale_global, slab_scale, error_scale);
  fitted += EV * beta_ev;
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
// RHS prior model
  target += std_normal_lpdf(z);
  target += std_normal_lpdf(aux1_local);
  target += inv_gamma_lpdf(aux2_local | 0.5, 0.5); // .5 * nu_local, .5 * nu_local, nu_local = 1
  target += std_normal_lpdf(aux1_global);
  target += inv_gamma_lpdf(aux2_global | 0.5, 0.5); // .5 * nu_local, .5 * nu_global, both = 1
  target += inv_gamma_lpdf(caux | 0.5*slab_df, 0.5*slab_df);
}

generated quantities {
  vector[n] esf;
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
      esf[i] = EV[i] * beta_ev;
#include parts/gen_quants_expression_in_loop.stan      
  }
}

