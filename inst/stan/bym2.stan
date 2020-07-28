functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=0> n_edges;
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  real<lower=0> scaling_factor; // scales the spatial component
}

transformed data {
#include parts/trans_data.stan
}

parameters {
  vector[n] phi;
  vector[n] theta;
  real<lower=0> sigma_re;
  real logit_rho;
#include parts/params.stan  
}

transformed parameters {
  real<lower=0, upper=1> rho;
  vector[n] convolved_re;
#include parts/trans_params_declaration.stan
  rho = inv_logit(logit_rho);
  convolved_re = sigma_re * (sqrt(rho / scaling_factor) * phi + sqrt(1 - rho) * theta);
  f += convolved_re;
#include parts/trans_params_expression.stan
}

model {
// BYM2 model
  phi ~ icar_normal(n, node1, node2);
  theta ~ std_normal();
  logit_rho ~ std_normal();
  sigma_re ~ std_normal();
#include parts/model.stan
}

generated quantities {
  vector[n] ssre; // scaled phi
  vector[n] sure; // scaled theta
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
    ssre[i] = sigma_re * sqrt(rho / scaling_factor) * phi[i];
    sure[i] = sigma_re * sqrt(1 - rho) * theta[i];
#include parts/gen_quants_expression_in_loop.stan      
  }
}

