// parameter models
  target += normal_lpdf(intercept | alpha_prior[1], alpha_prior[2]);
  if (dx_all) target += normal_lpdf(append_row(gamma, beta) | beta_prior[1], beta_prior[2]);  
  if (has_sigma) target += student_t_lpdf(sigma | sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_student) target += gamma_lpdf(nu[1] | t_nu_prior[1], t_nu_prior[2]);  
// data models (observational uncertainty)
  if (dx_me) {
    if (spatial_me) {
      for (j in 1:dx_me) {
        target += normal_lpdf(x_me[j] | x_true[j], sigma_me[j]);	       
	target += auto_normal_lpdf(x_true[j] |
				   rep_vector(mu_x_true[j], n),
				   sigma_x_true[j],
				   car_rho_x_true[j],
				   Ax_w, Ax_v, Ax_u,
				   Cidx,
				   Delta_inv, log_det_Delta_inv,
				   lambda, n, WCAR);	    
      }
    } else {
      for (j in 1:dx_me) {
	target += normal_lpdf(x_me[j] | x_true[j], sigma_me[j]);
	target += student_t_lpdf(x_true[j] | nu_x_true[j], mu_x_true[j], sigma_x_true[j]);
	target += gamma_lpdf(nu_x_true[j] | 3, 0.2);
      }
    }
      target += normal_lpdf(mu_x_true | prior_mux_true_location, prior_mux_true_scale);
      target += student_t_lpdf(sigma_x_true | 10, 0, prior_sigmax_true_scale);
}
// partial pooling of observations across all groups/geographies (varying intercept)
  if (has_re) {
    target += normal_lpdf(alpha_re | 0, alpha_tau[has_re]);
    target += student_t_lpdf(alpha_tau[has_re] | alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
  }
// process model (likelihood)
  if (!prior_only) {
    if (is_student)  target += student_t_lpdf(y | nu[1], fitted, sigma[has_sigma]);
    if (is_gaussian) target += normal_lpdf(y | fitted, sigma[has_sigma]);
    if (is_poisson)  target += poisson_lpmf(y_int | fitted);
    if (is_binomial) target += binomial_lpmf(y_int | trials, fitted);
}
  
