// parameter models
  array[dx_me] vector[n] x_true_transform;
  x_true_transform = x_true;
  target += normal_lpdf(intercept | prior_alpha[1], prior_alpha[2]);
  if (dx_all) target += normal_lpdf(append_row(gamma, beta) | prior_beta_location, prior_beta_scale);
  if (has_sigma) target += student_t_lpdf(sigma | prior_sigma[1], prior_sigma[2], prior_sigma[3]) - student_t_lcdf(0 | prior_sigma[1], prior_sigma[2], prior_sigma[3]);
  if (is_student) target += gamma_lpdf(nu[1] | prior_t_nu[1], prior_t_nu[2]);

// data models (observational uncertainty)
  if (dx_me) {
    if (spatial_me) {
      for (j in 1:dx_me) {
        vector[n] mu_x_true_tmp = rep_vector(mu_x_true[j], n);
        target += normal_lpdf(x_me[j] | x_true[j], sigma_me[j]);
	if (use_logit[j] > 0) x_true_transform[j] = logit(x_true[j]);
	target += auto_normal_lpdf(x_true_transform[j] |
				   mu_x_true_tmp,
				   sigma_x_true[j],
				   car_rho_x_true[j],
				   A_w,
				   A_v,
				   A_u,
				   Delta_inv,
				   log_det_Delta_inv,
				   lambda,
				   n,
				   WCAR);
      }
    } else {
      for (j in 1:dx_me) {
	target += normal_lpdf(x_me[j] | x_true[j], sigma_me[j]);
	if (use_logit[j] > 0) x_true_transform[j] = logit(x_true[j]);
	target += student_t_lpdf(x_true_transform[j] | nu_x_true[j], mu_x_true[j], sigma_x_true[j]);
	target += gamma_lpdf(nu_x_true[j] | prior_nux_true_alpha[j], prior_nux_true_beta[j]);
      }
    }
      target += normal_lpdf(mu_x_true | prior_mux_true_location, prior_mux_true_scale);
      target += student_t_lpdf(sigma_x_true | prior_sigmax_true_df, prior_sigmax_true_location, prior_sigmax_true_scale) - dx_me * student_t_lcdf(0 | prior_sigmax_true_df, prior_sigmax_true_location, prior_sigmax_true_scale);
}

// partial pooling of observations across all groups/geographies (varying intercept)
if (has_re) {
  target += normal_lpdf(alpha_re | 0, alpha_tau[has_re]);
  target += student_t_lpdf(alpha_tau[has_re] | prior_alpha_tau[1], prior_alpha_tau[2], prior_alpha_tau[3]) - student_t_lcdf(alpha_tau[has_re] | prior_alpha_tau[1], prior_alpha_tau[2], prior_alpha_tau[3]);
 }

// process model (likelihood)
if (!prior_only) {
  if (is_student)  target += student_t_lpdf(y[y_obs_idx] | nu[1], fitted[y_obs_idx], sigma[has_sigma]);
  if (is_gaussian) target += normal_lpdf(y[y_obs_idx] | fitted[y_obs_idx], sigma[has_sigma]);
  if (is_poisson) {
    target += poisson_lpmf(y_int[y_obs_idx] | fitted[y_obs_idx]);
    if (censor_point > 0) target += poisson_lcdf(censor_point | fitted[y_mis_idx]);
  }
  if (is_binomial) {
    target += binomial_lpmf(y_int[y_obs_idx] | trials[y_obs_idx], fitted[y_obs_idx]);
  }   
 }

  // ICAR
  if (type) {
    if (has_theta) {
      target += std_normal_lpdf(theta_tilde);
      if (type == 2) target += std_normal_lpdf(theta_scale[1]);
      // implicit uniform prior on rho:   if (type == 3) rho[1] ~ beta(1, 1);
      // no -lcdf(theta_scale[1]) because the ICAR model is improper anyways.
    }
    target += std_normal_lpdf(spatial_scale[1]);
    phi_tilde ~ icar_normal(spatial_scale[1], node1, node2, k, group_size, group_idx, has_theta);
    if (m) target += normal_lpdf(alpha_phi | 0, prior_alpha[2]);
  }

  // ESF
  if (dev) {
    target += std_normal_lpdf(z);
    target += std_normal_lpdf(aux1_local) - std_normal_lcdf(0);
    target += inv_gamma_lpdf(aux2_local | 0.5, 0.5); // .5 * nu_local, .5 * nu_local, nu_local = 1
    target += std_normal_lpdf(aux1_global[1]) - std_normal_lcdf(0);
    target += inv_gamma_lpdf(aux2_global[1] | 0.5, 0.5); // .5 * nu_local, .5 * nu_global, both = 1
    target += inv_gamma_lpdf(caux[1] | 0.5*slab_df, 0.5*slab_df);
  }

  // CAR
if (car) {
    
  target += student_t_lpdf(car_scale[1] | prior_sigma[1], prior_sigma[2], prior_sigma[3]) - student_t_lcdf(0 | prior_sigma[1], prior_sigma[2], prior_sigma[3]);
  
  if (is_auto_gaussian && prior_only == 0) {
    target += auto_normal_lpdf(y |
			       fitted,
			       car_scale[1],
			       car_rho[1],
			       A_w,
			       A_v,
			       A_u,
			       Delta_inv,
			       log_det_Delta_inv,
			       lambda,
			       n,
			       WCAR);
  }
  
  if (!is_auto_gaussian) {

    if (ZMP) {
      // zero-mean constrained/parameterized CAR model (hierarchical)
      	target += auto_normal_lpdf(log_lambda |
				   zero_vec,
				   1, // scale
				   car_rho[1],
				   A_w,
				   A_v,
				   A_u,
				   Delta_inv,
				   log_det_Delta_inv,
				   lambda,
				   n,
				   WCAR);
      
    } else {

      // classic (modeled mean) hierarchical CAR model 
      	target += auto_normal_lpdf(log_lambda |
				   log_lambda_mu,
				   car_scale[1],
				   car_rho[1],
				   A_w,
				   A_v,
				   A_u,
				   Delta_inv,
				   log_det_Delta_inv,
				   lambda,
				   n,
				   WCAR);
    }
  }
    
 }

 // SAR
if (sar > 0) {
  target += student_t_lpdf(sar_scale[1] | prior_sigma[1], prior_sigma[2], prior_sigma[3]) - student_t_lcdf(0 | prior_sigma[1], prior_sigma[2], prior_sigma[3]);
  if (is_auto_gaussian && prior_only == 0) {
    target += sar_normal_lpdf(y |
			      fitted,
			      sar_scale[1],
			      sar_rho[1],
			      W_w,
			      W_v,
			      W_u,
			      eigenvalues_w,
			      n,
			      sar);
  }
  
  if (!is_auto_gaussian) {

    if (ZMP) {
    
      target += sar_normal_lpdf(log_lambda |
				zero_vec,
				1, //scale
				sar_rho[1],
				W_w,
				W_v,
				W_u,
				eigenvalues_w,
				n,
				sar);
      
    } else {
    
    target += sar_normal_lpdf(log_lambda |
			      log_lambda_mu,
			      sar_scale[1],
			      sar_rho[1],
			      W_w,
			      W_v,
			      W_u,
			      eigenvalues_w,
			      n,
			      sar);			      
    }
  }
 }

