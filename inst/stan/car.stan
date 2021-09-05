functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=0, upper=1> invert;
  real car_rho_lims[2];
  matrix<lower=0>[invert ? n : 1, invert ? n : 1] C;  
}

transformed data {
#include parts/trans_data.stan
}

parameters {
  vector[is_auto_gaussian ? 0 : n] log_lambda;
  real<lower=0> car_scale;
  real<lower=car_rho_lims[1], upper=car_rho_lims[2]> car_rho;   
#include parts/params.stan
}

transformed parameters {
  // declaration
  matrix[n, dx_all] x_all;
  vector[n] log_lambda_mu;
  vector[n] fitted;  
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me) for (j in 1:dx_me) x_all[ ,x_me_idx[j]] = x_true[j];
  // car
  log_lambda_mu = rep_vector(intercept, n);
  if (has_re) {
    for (i in 1:n) {
      log_lambda_mu[i] += alpha_re[id[i]];
    }
  }  
  if (dwx) {
    if (has_me) {
      for (i in 1:dwx) {
	log_lambda_mu += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];
      }
    } else {
      log_lambda_mu += WX * gamma;
    }
  } 
  if (dx_all) log_lambda_mu += x_all * beta;
  if (is_auto_gaussian) {
    fitted = offset + log_lambda_mu;
      } else {
    fitted = offset + log_lambda;
  }
  if (is_binomial) fitted = inv_logit(fitted);
  if (is_poisson) fitted = exp(fitted);
}

model {
#include parts/model.stan
  target += student_t_lpdf(car_scale | sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_auto_gaussian * !prior_only) {
      target += auto_normal_lpdf(y |
				 fitted, car_scale, car_rho,
				 Ax_w, Ax_v, Ax_u,
				 Cidx,
				 Delta_inv, log_det_Delta_inv,
				 lambda, n, WCAR);
  }
  if (!is_auto_gaussian) {
    target += auto_normal_lpdf(log_lambda |
				log_lambda_mu, car_scale, car_rho,
				Ax_w, Ax_v, Ax_u,
				Cidx,
				Delta_inv, log_det_Delta_inv,
			       lambda, n, WCAR);
    }
}

generated quantities {
  matrix[invert ? n : 0, invert ? n : 0] S;
  vector[is_auto_gaussian ? 0 : n] phi;
#include parts/gen_quants_declaration.stan
  if (!is_auto_gaussian) phi = log_lambda - log_lambda_mu;  
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan
  }
  if (invert * is_auto_gaussian) {
    S = car_scale^2 * diag_post_multiply(inverse(diag_matrix(rep_vector(1, n)) - car_rho * C), 1.0 ./ Delta_inv);      
    log_lik[1] = multi_normal_lpdf(y | fitted, S);
  } 
}

