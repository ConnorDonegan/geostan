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
  fitted = input_offset;
#include parts/trans_params_expression_icar.stan
#include parts/trans_params_expression_esf.stan
#include parts/trans_params_expression_auto-model.stan
  
  if (as_olm) {

    fitted += intercept;
    
    if (has_re) fitted += alpha_re[id];
    
    if (use_qr) {

      fitted += Q_ast * coefs_qr;
      coefs = R_ast_inverse * coefs_qr;
      beta = coefs[1:dx_all];
      if (dwx) gamma = coefs[(dx_all+1):d_qr];	       
      
    } else {
      
      if (dx_all) {
	beta = coefs_qr[1:dx_all];	
	fitted += x_all * beta;
      }
      
      if (dwx) {
	gamma = coefs_qr[(dx_all+1):d_qr];	
	for (i in 1:dwx) fitted += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];
      }
      
    }
    
  }
  
  if (is_binomial) fitted = inv_logit(fitted);
  if (is_poisson) fitted = exp(fitted);
  
}

model {
#include parts/model.stan
}


