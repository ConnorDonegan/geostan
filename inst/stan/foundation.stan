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
#include parts/trans_params_expression_auto-model.stan
  if (car == 0 && sar == 0) {
    
    if (has_re) fitted += alpha_re[id];
    
    if (use_qr) {      
      fitted += Q_ast[ , 1:dx_obs] * beta_qr;
      beta = R_ast_inverse[1:dx_obs, 1:dx_obs]  * beta_qr;
      if (dwx) {
	fitted += Q_ast[ , (dx_obs+1):d_qr] * gamma_qr;
	gamma = R_ast_inverse[(dx_obs+1):d_qr, (dx_obs+1):d_qr] * gamma_qr;
      }      
    } else {      
      if (dwx) {
	for (i in 1:dwx) {
	  fitted += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma_qr[i];
	}
	gamma = gamma_qr;
      }      
      if (dx_all) {
	fitted += x_all * beta_qr;
	beta = beta_qr;
      }      
    }    
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

