functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  matrix<lower=0, upper=1>[n, n] C; // adjacency matrix
  int<lower=0> dc_nonzero; // number of non-zero elements in c
  vector[n] D_diag; // diagonal of D, e.g., number of neighbors
  int<lower=0, upper=1> invert;
}

transformed data {
  vector[n] mean_zero = rep_vector(0, n);
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  vector[n] invsqrtD;
  vector[dc_nonzero] car_w = csr_extract_w(C'); 
  int car_v[dc_nonzero] = csr_extract_v(C');    
  int car_u[n + 1] = csr_extract_u(C');         
#include parts/trans_data.stan
  for (i in 1:n) invsqrtD[i] = 1 / sqrt(D_diag[i]);
  lambda = eigenvalues_sym(quad_form_diag(C, invsqrtD));
}

parameters {
  vector[is_auto_gaussian ? 0 : n] phi;
  real<lower = 0> car_precision;
  real<lower=1/min(lambda), upper=1/max(lambda)> car_alpha;   
#include parts/params.stan
}

transformed parameters {
  real<lower=0> car_scale = 1 / sqrt(car_precision);
#include parts/trans_params_declaration.stan
  if (!is_auto_gaussian) f += phi;
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
  car_scale ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_auto_gaussian * !prior_only) y ~ car_normal(f, car_precision, car_alpha, car_w, car_v, car_u, D_diag, lambda, n);
  if (!is_auto_gaussian * !prior_only) phi ~ car_normal(mean_zero, car_precision, car_alpha, car_w, car_v, car_u, D_diag, lambda, n);
}

generated quantities {
  matrix[n, n] S;
  vector[n] trend;
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan
  }
  if (is_auto_gaussian) {
  trend = car_alpha * C * (y - f) ./ D_diag;
  fitted = f;
  residual = y - f - trend;
 }
  
  if (invert * is_auto_gaussian) {
    S = car_scale^2 * inverse(diag_matrix(D_diag) - car_alpha * C);      
    yrep = multi_normal_rng(f, S);    
    log_lik[1] = multi_normal_lpdf(y | f, S);
  }
}

