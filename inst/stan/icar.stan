functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=0> n_edges; 
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges]; 
  real<lower=0> phi_scale_prior; 
}

transformed data {
#include parts/trans_data.stan
}

parameters {
  vector[n] phi_tilde;
  real<lower=0> phi_scale;
#include parts/params.stan
}

transformed parameters {
#include parts/trans_params_declaration.stan
  f += phi_tilde * phi_scale;
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
// IAR model
  phi_tilde ~ icar_normal(n, node1, node2);
  phi_scale ~ normal(0, phi_scale_prior);
 }

generated quantities {
  vector[n] phi;
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
      phi[i] = phi_tilde[i] * phi_scale;
#include parts/gen_quants_expression_in_loop.stan      
  }
}

