// parameter models
  target += normal_lpdf(intercept | alpha_prior[1], alpha_prior[2]);
  if (dx_all) target += normal_lpdf(append_row(gamma, beta) | beta_prior[1], beta_prior[2]);  
  if (has_sigma) target += student_t_lpdf(sigma | sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_student) target += gamma_lpdf(nu[1] | t_nu_prior[1], t_nu_prior[2]);  
// data models (observational uncertainty)
  // unbounded variables
  if (dx_me_unbounded) {
    if (spatial_me) {
      for (j in 1:dx_me_unbounded) {
        target += normal_lpdf(x_me_unbounded[j] | x_true_unbounded[j], sigma_me_unbounded[j]);	       
	if (WCAR) {	  
	  target += wcar_normal_lpdf(x_true_unbounded[j] |
				      rep_vector(mu_x_true_unbounded[j], n),
				      sigma_x_true_unbounded[j],
				      car_rho_x_true_unbounded[j],
				      Ax_w, Ax_v, Ax_u,
				      Delta_inv, log_det_Delta_inv,
				      lambda, n);
	} else {
	    target += car_normal_lpdf(x_true_unbounded[j] |
		       rep_vector(mu_x_true_unbounded[j], n),
		       sigma_x_true_unbounded[j],
		       car_rho_x_true_unbounded[j],
		       Ax_w, Ax_v, Ax_u,
		       Cidx,
		       Delta_inv, log_det_Delta_inv,
		       lambda, n);	    
	}
      }
      } else {
      for (j in 1:dx_me_unbounded) {
	target += normal_lpdf(x_me_unbounded[j] | x_true_unbounded[j], sigma_me_unbounded[j]);
	target += student_t_lpdf(x_true_unbounded[j] | nu_x_true_unbounded[j], mu_x_true_unbounded[j], sigma_x_true_unbounded[j]);
	target += gamma_lpdf(nu_x_true_unbounded[j] | 3, 0.2);
    }
  }
      target += normal_lpdf(mu_x_true_unbounded | prior_mean_x_true_unbounded, 2 * prior_scale_x_true_unbounded);
      target += student_t_lpdf(sigma_x_true_unbounded | 10, 0, 2 * prior_scale_x_true_unbounded);
}
  // bounded variables
  if (dx_me_bounded) {
    if (spatial_me) {
      for (j in 1:dx_me_bounded) {
	target += normal_lpdf(x_me_bounded[j] | x_true_bounded[j], sigma_me_bounded[j]);
	if (WCAR) {
	   target += wcar_normal_lpdf(x_true_bounded[j] |
				      rep_vector(mu_x_true_bounded[j], n),
				      sigma_x_true_bounded[j],
				      car_rho_x_true_bounded[j],
				      Ax_w, Ax_v, Ax_u,
				      Delta_inv, log_det_Delta_inv,
				      lambda, n);
	} else {
	  target += car_normal_lpdf(x_true_bounded[j] |
				    rep_vector(mu_x_true_bounded[j], n),
				    sigma_x_true_bounded[j],
				    car_rho_x_true_bounded[j],
				    Ax_w, Ax_v, Ax_u,
				    Cidx,
				    Delta_inv, log_det_Delta_inv,
				    lambda, n);					 
	}	
        }
      } else {
      for (j in 1:dx_me_bounded) {
	target += normal_lpdf(x_me_bounded[j] | x_true_bounded[j], sigma_me_bounded[j]);		
	target += student_t_lpdf(x_true_bounded[j] | nu_x_true_bounded[j], mu_x_true_bounded[j], sigma_x_true_bounded[j]);
	target += gamma_lpdf(nu_x_true_bounded[j] | 3, 0.2);
    }
  }
    target += normal_lpdf(mu_x_true_bounded | prior_mean_x_true_bounded, 2 * prior_scale_x_true_bounded);
    target += student_t_lpdf(sigma_x_true_bounded | 10, 0, 2 * prior_scale_x_true_bounded);
}
// partial pooling of observations across all groups/geographies (varying intercept)
  if (has_re) {
    target += student_t_lpdf(alpha_tau[has_re] | alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    target += std_normal_lpdf(alpha_re_tilde);
  }
// process model (likelihood)
  if (!prior_only) {
    if (is_student) {	
      target += student_t_lpdf(y | nu[1], fitted, sigma[has_sigma]);
      }   
    if (is_gaussian) target += normal_lpdf(y | fitted, sigma[has_sigma]);
    if (is_poisson) target += poisson_lpmf(y_int | fitted);
    if (is_binomial) target += binomial_lpmf(y_int | trials, fitted);
}
  
