// parameter models 
  intercept ~ normal(alpha_prior[1], alpha_prior[2]);
  if (dx_all) append_row(gamma, beta) ~ normal(beta_prior[1], beta_prior[2]);
  if (has_sigma) sigma ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_student) nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);  
// data models (observational uncertainty)
  // unbounded variables
  if (dx_me_unbounded) {
    if (spatial_me) {
      for (j in 1:dx_me_unbounded) {
        x_me_unbounded[j] ~ normal(x_true_unbounded[j], sigma_me_unbounded[j]);
	x_true_unbounded[j] ~ car_normal(rep_vector(mu_x_true_unbounded[j], n), sigma_x_true_unbounded[j], car_rho_x_true_unbounded[j], ImC, ImC_v, ImC_u, Cidx, M_inv, lambda, n);
        }
      } else {
      for (j in 1:dx_me_unbounded) {
      	  x_me_unbounded[j] ~ normal(x_true_unbounded[j], sigma_me_unbounded[j]);
	  x_true_unbounded[j] ~ student_t(nu_x_true_unbounded[j], mu_x_true_unbounded[j], sigma_x_true_unbounded[j]);
	  nu_x_true_unbounded[j] ~ gamma(3, 0.2);
    }
  }
   mu_x_true_unbounded ~ normal(prior_mean_x_true_unbounded, 2 * prior_scale_x_true_unbounded);
   sigma_x_true_unbounded ~ student_t(10, 0, 2 * prior_scale_x_true_unbounded);
}
  // bounded variables
  if (dx_me_bounded) {
    if (spatial_me) {
      for (j in 1:dx_me_bounded) {
	x_me_bounded[j] ~ normal(x_true_bounded[j], sigma_me_bounded[j]);	
        x_true_bounded[j] ~ car_normal(rep_vector(mu_x_true_bounded[j], n), sigma_x_true_bounded[j], car_rho_x_true_bounded[j], ImC, ImC_v, ImC_u, Cidx, M_inv, lambda, n);
        }
      } else {
      for (j in 1:dx_me_bounded) {
	x_me_bounded[j] ~ normal(x_true_bounded[j], sigma_me_bounded[j]);		
	x_true_bounded[j] ~ student_t(nu_x_true_bounded[j], mu_x_true_bounded[j], sigma_x_true_bounded[j]);
	nu_x_true_bounded[j] ~ gamma(3, 0.2);
    }
  }
   mu_x_true_bounded ~ normal(prior_mean_x_true_bounded, 2 * prior_scale_x_true_bounded);
   sigma_x_true_bounded ~ student_t(10, 0, 2 * prior_scale_x_true_bounded);
}
// partial pooling of observations across all groups/geographies (varying intercept)
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
// process model (likelihood)
  if (!prior_only) {
    if (is_student) {	
      y ~ student_t(nu[1], fitted, sigma[has_sigma]);
      }   
    if (is_gaussian)  y ~ normal(fitted, sigma[has_sigma]);
    if (is_poisson) y_int ~ poisson(fitted);
    if (is_binomial) y_int ~ binomial(trials, fitted);
}
  
