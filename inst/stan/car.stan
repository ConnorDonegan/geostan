functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=0, upper=1> invert;
}

transformed data {
#include parts/trans_data.stan
}

parameters {
  vector[is_auto_gaussian ? 0 : n] phi;
  real<lower=0> car_scale;
  real<lower=1/min(lambda), upper=1/max(lambda)> car_rho;   
#include parts/params.stan
}

transformed parameters {
  // declaration
  matrix[n, dx_all] x_all;
  vector[n] f;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me_unbounded) for (j in 1:dx_me_unbounded) x_all[ ,x_me_unbounded_idx[j]] = x_true_unbounded[j];
  if (dx_me_bounded) for (j in 1:dx_me_bounded) x_all[,x_me_bounded_idx[j]] = x_true_bounded[j];
  f = offset; // add intercept below.
  // car
  if (is_auto_gaussian) f += intercept;
  if (!is_auto_gaussian) f += phi; // phi ~ Gau(intercpet, (I - rho C)^-1 M)
  // expression
    if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
   }
  }  
  if (dwx) {
   if (has_me) {
      for (i in 1:dwx) {
     f += csr_matrix_times_vector(n, n, w, v, u, x_all[,wx_idx[i]]) * gamma[i];
     }
   } else {
      f += WX * gamma;
      }
  } 
  if (dx_all) f += x_all * beta;
  if (is_binomial) f = inv_logit(f);

  //#include parts/trans_params_declaration.stan
  //  if (!is_auto_gaussian) f += phi;
  //#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
  car_scale ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_auto_gaussian * !prior_only) y ~ car_normal(f, car_scale, car_rho, ImC, ImC_v, ImC_u, Cidx, M_inv, lambda, n);
  if (!is_auto_gaussian) phi ~ car_normal(rep_vector(intercept, n), car_scale, car_rho, ImC, ImC_v, ImC_u, Cidx, M_inv, lambda, n);
}

generated quantities {
  matrix[n, n] S;
  vector[n] trend;
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan
  }
  if (is_auto_gaussian) {
  trend = car_rho * C * (y - f);
  fitted = f;
  residual = y - f - trend;
 }
  
  if (invert * is_auto_gaussian) {
    S = car_scale^2 * diag_post_multiply(inverse(diag_matrix(rep_vector(1, n)) - car_rho * C), M_diag);      
    yrep = multi_normal_rng(f, S);    
    log_lik[1] = multi_normal_lpdf(y | f, S);
  }
}

