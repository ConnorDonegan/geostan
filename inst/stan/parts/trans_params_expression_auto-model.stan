if (car > 0 || sar > 0) {
  if (is_auto_gaussian == 0) {
    // if auto-Normal model: treat as ordinary linear model (OLMish)
    
    if (ZMP) {
      // Hierarchical model with SAR or CAR model centered on zero
      // 'zero-mean parameterized' (ZMP) spatial model

        // log_lambda ~ MVN(0, 1 * Sigma);
      if (car > 0) fitted += log_lambda * car_scale[1];
      if (sar > 0) fitted += log_lambda * sar_scale[1];
        //log_lambda_mu = zero_vec;
	
    }
  
    if (ZMP == 0) {
      // MMP: modeled-mean parameterization (MMP) of CAR or SAR model for (n-vector) rate parameter     
      // Hierarchical CAR or SAR model ONLY  (must be Poisson or binomial model for Y)

      // log_lambda_mu should only be declared if NCP==0 and (SAR||CAR) == 1
      log_lambda_mu = rep_vector(intercept, n);
     
      if (has_re) log_lambda_mu += alpha_re[id];
     
      if (use_qr) {
	
	log_lambda_mu += Q_ast * coefs_qr;
	coefs = R_ast_inverse * coefs_qr;
	beta = coefs[1:dx_all];
	if (dwx) gamma = coefs[(dx_all+1):d_qr];	       	
	
      } else {

	if (dx_all) {
	  beta = coefs_qr[1:dx_all];		  
	  log_lambda_mu += x_all * beta;
	}
	
	if (dwx) {
	  gamma = coefs_qr[(dx_all+1):d_qr];		  
	  for (i in 1:dwx) log_lambda_mu += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];	  
	}	
	
      }
      
      fitted += log_lambda;
       // log_lambda ~ MVN(log_lambda_mu, tau^2 * Sigma);      
    }
    
  }
 }

