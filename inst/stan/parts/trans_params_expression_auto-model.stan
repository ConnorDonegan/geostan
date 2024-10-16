  if (car > 0 || sar > 0) {
    
    log_lambda_mu = rep_vector(intercept, n);
    
    if (has_re) log_lambda_mu += alpha_re[id];

    if (use_qr) {
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
      
    if (is_auto_gaussian) {
      fitted = input_offset + log_lambda_mu;
    } else {
      fitted = input_offset + log_lambda;
    }
    
  }
