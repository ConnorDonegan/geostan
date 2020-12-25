// spatial ME model parameter declaration
  // unbounded
  vector[dev] beta_ev_me_unbounded[spatial_me ? dx_me_unbounded : 0];
  vector[n] fitted_me_unbounded[spatial_me ? dx_me_unbounded : 0];
  // bounded
  vector[dev] beta_ev_me_bounded[spatial_me ? dx_me_bounded : 0];
  vector[n] fitted_me_bounded[spatial_me ? dx_me_bounded : 0];
// parameter models 
  intercept ~ normal(alpha_prior[1], alpha_prior[2]);
  if (dx_all) append_row(gamma, beta) ~ normal(beta_prior[1], beta_prior[2]);
  if (has_sigma) sigma ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]);
  if (is_student) nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);  
// data models (observational error)
  if (model_offset) {
    offset_obs ~ normal(offset_est, offset_me);
    offset_est ~ student_t(nu_offset[1], mu_offset[1], sigma_offset[1]);
    nu_offset[1] ~ gamma(3, 0.2);
  }
  if (dx_me_unbounded) {
    if (spatial_me) {
      for (j in 1:dx_me_unbounded) {
        beta_ev_me_unbounded[j] = rhs_prior(dev, z_ev_me_unbounded[j],
	                                    aux1_global_me_unbounded[j],
					    aux2_global_me_unbounded[j],
					    aux1_local_me_unbounded[j],
					    aux2_local_me_unbounded[j],
					    caux_me_unbounded[j],
					    scale_global_me_unbounded[j],
					    slab_scale_me_unbounded[j],
					    sigma_x_true_unbounded[j]);
        fitted_me_unbounded[j] = mu_x_true_unbounded[j] + EV * beta_ev_me_unbounded[j];
        z_ev_me_unbounded[j] ~ std_normal();
        aux1_local_me_unbounded[j] ~ std_normal();
        aux2_local_me_unbounded[j] ~ inv_gamma(0.5, 0.5);
        caux_me_unbounded ~ inv_gamma(0.5 * slab_df_me_unbounded[j], 0.5 * slab_df_me_unbounded[j]);		
        x_me_unbounded[j] ~ normal(x_true_unbounded[j], sigma_me_unbounded[j]);	
        x_true_unbounded[j] ~ student_t(nu_x_true_unbounded[j], fitted_me_unbounded[j], sigma_x_true_unbounded[j]);	
        }
      aux1_global_me_unbounded ~ std_normal();
      aux2_global_me_unbounded ~ inv_gamma(0.5, 0.5);
      } else {
      for (j in 1:dx_me_unbounded) {
      	  x_me_unbounded[j] ~ normal(x_true_unbounded[j], sigma_me_unbounded[j]);
	  x_true_unbounded[j] ~ student_t(nu_x_true_unbounded[j], mu_x_true_unbounded[j], sigma_x_true_unbounded[j]);
    }
  }
   nu_x_true_unbounded ~ gamma(3, 0.2);   
}
  if (dx_me_bounded) {
    if (spatial_me) {
      for (j in 1:dx_me_bounded) {
        beta_ev_me_bounded[j] = rhs_prior(dev, z_ev_me_bounded[j],
	                                    aux1_global_me_bounded[j],
					    aux2_global_me_bounded[j],
					    aux1_local_me_bounded[j],
					    aux2_local_me_bounded[j],
					    caux_me_bounded[j],
					    scale_global_me_bounded[j],
					    slab_scale_me_bounded[j],
					    sigma_x_true_bounded[j]);
        fitted_me_bounded[j] = mu_x_true_bounded[j] + EV * beta_ev_me_bounded[j];
        z_ev_me_bounded[j] ~ std_normal();
        aux1_local_me_bounded[j] ~ std_normal();
        aux2_local_me_bounded[j] ~ inv_gamma(0.5, 0.5);
        caux_me_bounded ~ inv_gamma(0.5 * slab_df_me_bounded[j], 0.5 * slab_df_me_bounded[j]);
        x_me_bounded[j] ~ normal(x_true_bounded[j], sigma_me_bounded[j]);	
        x_true_bounded[j] ~ student_t(nu_x_true_bounded[j], fitted_me_bounded[j], sigma_x_true_bounded[j]);
        }
      aux1_global_me_bounded ~ std_normal();
      aux2_global_me_bounded ~ inv_gamma(0.5, 0.5);
      } else {
      for (j in 1:dx_me_bounded) {
      	  x_me_bounded[j] ~ normal(x_true_bounded[j], sigma_me_bounded[j]);
	  x_true_bounded[j] ~ student_t(nu_x_true_bounded[j], mu_x_true_bounded[j], sigma_x_true_bounded[j]);
    }
  }
   nu_x_true_bounded ~ gamma(3, 0.2);   
}
// partial pooling of observations across groups/geographies (exchangeable "random effects")
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
// process model (likelihood of the data)
  if (!prior_only) {
    if (is_student) {	
      y ~ student_t(nu[1], f, sigma[has_sigma]);
      }   
    if (is_gaussian)  y ~ normal(f, sigma[has_sigma]);
    if (is_poisson) y_int ~ poisson_log(f);
    if (is_binomial) y_int ~ binomial(trials, f);
  }
  
