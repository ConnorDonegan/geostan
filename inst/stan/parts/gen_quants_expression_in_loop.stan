   if (is_student) {
      fitted[i] = f[i];
      residual[i] = y[i] - fitted[i];
      log_lik[i] = student_t_lpdf(y[i] | nu[1], fitted[i], sigma[has_sigma]);
      yrep[i] = student_t_rng(nu[1], fitted[i], sigma[has_sigma]);     
   }
   if (is_gaussian) {
         fitted[i] = f[i];
	 residual[i] = y[i] - fitted[i];
      log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma[has_sigma]);
      yrep[i] = normal_rng(fitted[i], sigma[has_sigma]);     
   }
   if (is_poisson) {
       fitted[i] = exp(f[i] - offset[i]);
       residual[i] = ( y[i] / exp(offset[i]) ) - fitted[i];
       log_lik[i] = poisson_log_lpmf(y_int[i] | f[i]);
       if (f[i] > 20) {
       	  print("f[i] too large (>20) for poisson_log_rng");
       	  yrep[i] = -1;
       } else {
       	 yrep[i] = poisson_log_rng(f[i]);
   }
  }
  if (is_binomial) {
    fitted[i] = (f[i] * trials[i]) / trials[i];
    residual[i] = (y[i] / trials[i]) - fitted[i]; 
    yrep[i] = binomial_rng(trials[i], f[i]);
    log_lik[i] = binomial_lpmf(y_int[i] | trials[i], f[i]);
  }  



