functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  matrix<lower=0, upper=1>[n, n] C; // Binary connectivity matrix
  int C_n; // nrow(edges), number of edges or connected nodes
  int<lower=1, upper=n> C_sparse[C_n, 2]; // the edges, unique pairs of connected nodes (i < j)
  vector<lower=0>[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  real<lower=0> phi_tau_prior[2];
}

transformed data {
  vector[n] lambda;       // eigenvalues of invsqrtD * C * invsqrtD
  vector[n] invsqrtD;    
#include parts/trans_data.stan
  for (i in 1:n) invsqrtD[i] = 1 / sqrt(D_sparse[i]);
  lambda = eigenvalues_sym(quad_form(C, diag_matrix(invsqrtD)));
}

parameters {
  vector[n] phi;
  real<lower = 0> phi_tau;
  real<lower = 0, upper = 1> phi_alpha;   // implicit uniform(0, 1) prior
#include parts/params.stan
}

transformed parameters {
#include parts/trans_params_declaration.stan
  f += phi;
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
// CAR model
  phi ~ sparse_car(phi_tau, phi_alpha, C_sparse, D_sparse, lambda, n, C_n);
  phi_tau ~ gamma(phi_tau_prior[1], phi_tau_prior[2]);
 }

generated quantities {
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan      
  }
}

