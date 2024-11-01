if (car || sar) {
  if (is_auto_gaussian == 0) {
    // if auto-Normal model: treat as ordinary linear model (OLMish)

    
    if (ZMP) {
      // Hierarchical model with SAR or CAR model centered on zero
      // 'zero-mean parameterized' (ZMP) spatial model

      // log_lambda ~ MVN(0, 1 * Sigma);
      if (car == 1) fitted += log_lambda * car_scale[1];
      if (sar == 1) fitted += log_lambda * sar_scale[1];
      //log_lambda_mu = zero_vec;
	
    }
  
    if (ZMP == 0) {
      // MMP: modeled-mean parameterization (MMP) of CAR or SAR model for (n-vector) rate parameter     
      // Hierarchical CAR or SAR model ONLY  (must be Poisson or binomial model for Y)

      // log_lambda_mu should only be declared if NCP==0 and (SAR||CAR) == 1
      log_lambda_mu = rep_vector(intercept, n);
     
      if (has_re) log_lambda_mu += alpha_re[id];
     
      if (use_qr) { //use qr
	
	log_lambda_mu += Q_ast[ , 1:dx_obs] * beta_qr;
	beta = R_ast_inverse[1:dx_obs, 1:dx_obs]  * beta_qr;
	
	if (dwx) {
	  log_lambda_mu += Q_ast[ , (dx_obs+1):d_qr] * gamma_qr;
	  gamma = R_ast_inverse[(dx_obs+1):d_qr, (dx_obs+1):d_qr] * gamma_qr;
	}
	
      } else {
	
	if (dwx) {
	  for (i in 1:dwx) log_lambda_mu += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma_qr[i];
	  gamma = gamma_qr;
	}
	
	if (dx_all) {
	  log_lambda_mu += x_all * beta_qr;
	  beta = beta_qr;
	}
	
      }
      
      fitted += log_lambda;
     // log_lambda ~ MVN(log_lambda_mu, tau^2 * Sigma);      
    }
    
  }
 }


/// current method
/* if (car || sar) { */
    
/*     log_lambda_mu = rep_vector(intercept, n); */
    
/*     if (has_re) log_lambda_mu += alpha_re[id]; */

/*     if (use_qr) { */
/*       log_lambda_mu += Q_ast[ , 1:dx_obs] * beta_qr; */
/*       beta = R_ast_inverse[1:dx_obs, 1:dx_obs]  * beta_qr; */
/*       if (dwx) { */
/* 	log_lambda_mu += Q_ast[ , (dx_obs+1):d_qr] * gamma_qr; */
/* 	gamma = R_ast_inverse[(dx_obs+1):d_qr, (dx_obs+1):d_qr] * gamma_qr; */
/*       } */
/*     } else {     */
/*       if (dwx) { */
/* 	for (i in 1:dwx) log_lambda_mu += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma_qr[i]; */
/* 	gamma = gamma_qr; */
/*       }     */
/*       if (dx_all) { */
/* 	log_lambda_mu += x_all * beta_qr; */
/* 	beta = beta_qr; */
/*       } */
/*     } */
      
/*     if (is_auto_gaussian) { */      
/*       // y ~ CAR or y ~ SAR */
/*       fitted += log_lambda_mu; */
      
/*     } else { */       
/*       fitted += log_lambda; */      
/*       } */
/*     } */
    
/*   } */

