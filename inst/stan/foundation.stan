functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
}

transformed data {
#include parts/trans_data.stan
}

parameters {
#include parts/params.stan
}

transformed parameters {
#include parts/trans_params_declaration.stan
#include parts/trans_params_expression_icar.stan
#include parts/trans_params_expression_esf.stan
//#include parts/trans_params_expression_car.stan
#include parts/trans_params_expression_auto-model.stan
  if (!car && !sar) {
    if (has_re) {
      for (i in 1:n) {
	fitted[i] += alpha_re[id[i]];
      }
    }  
    if (dwx) {
      for (i in 1:dwx) fitted += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];      
    } 
    if (dx_all) fitted += x_all * beta;
  }
  if (is_binomial) fitted = inv_logit(fitted);
  if (is_poisson) fitted = exp(fitted);  
}

model {
#include parts/model.stan
}

generated quantities {
#include parts/gen_quants_declaration.stan
#include parts/gen_quants_expression.stan
}

