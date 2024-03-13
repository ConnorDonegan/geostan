  if (is_auto_gaussian && car) {
    log_lik[1] = auto_normal_lpdf(y |
				  fitted, car_scale[1], car_rho[1],
				  Ax_w, Ax_v, Ax_u,
				  Cidx,
				  Delta_inv, log_det_Delta_inv,
				  lambda, n, WCAR);    
  }
  if (is_auto_gaussian && sar) {
    log_lik[1] = sar_normal_lpdf(y |
				 fitted,
				 sar_scale[1],
				 sar_rho[1],
				 ImW_w,
				 ImW_v,
				 ImW_u,
				 Widx,
				 eigenvalues_w,
				 n);
  }
  for (i in 1:n_obs) {
   if (is_student) {
      log_lik[i] = student_t_lpdf(y[y_obs_idx[i]] | nu[1], fitted[y_obs_idx[i]], sigma[has_sigma]); 
   }
   if (is_gaussian) {
      log_lik[i] = normal_lpdf(y[y_obs_idx[i]] | fitted[y_obs_idx[i]], sigma[has_sigma]);   
   }
   if (is_poisson) {
     log_lik[i] = poisson_lpmf(y_int[y_obs_idx[i]] | fitted[y_obs_idx[i]]);
  }
  if (is_binomial) {
    log_lik[i] = binomial_lpmf(y_int[y_obs_idx[i]] | trials[y_obs_idx[i]], fitted[y_obs_idx[i]]);
  }   
  }
