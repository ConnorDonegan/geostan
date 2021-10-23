  if (is_auto_gaussian) {
    log_lik[1] = auto_normal_lpdf(y |
				  fitted, car_scale[1], car_rho[1],
				  Ax_w, Ax_v, Ax_u,
				  Cidx,
				  Delta_inv, log_det_Delta_inv,
				  lambda, n, WCAR);    
  }
  for (i in 1:n) {
   if (is_student) {
      log_lik[i] = student_t_lpdf(y[i] | nu[1], fitted[i], sigma[has_sigma]); 
   }
   if (is_gaussian) {
      log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma[has_sigma]);   
   }
   if (is_poisson) {
       log_lik[i] = poisson_lpmf(y_int[i] | fitted[i]);
  }
  if (is_binomial) {
     log_lik[i] = binomial_lpmf(y_int[i] | trials[i], fitted[i]); 
  }   
  }
