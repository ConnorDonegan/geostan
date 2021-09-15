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
