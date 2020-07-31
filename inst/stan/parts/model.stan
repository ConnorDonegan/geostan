// parameter models (priors)
  intercept ~ normal(alpha_prior[1], alpha_prior[2]);
  if (dx_all) append_row(gamma, beta) ~ normal(beta_prior[1], beta_prior[2]);
  if (has_sigma) sigma ~ student_t(sigma_prior[1], sigma_prior[2], sigma_prior[3]);
// data models 
  if (model_offset) offset_obs ~ normal(offset_est, offset_me);
  if (dx_me_bounded) to_vector(x_me_bounded) ~ normal(to_vector(x_true_bounded), to_vector(sigma_me_bounded));
  if (dx_me_unbounded) to_vector(x_me_unbounded) ~ normal(to_vector(x_true_unbounded), to_vector(sigma_me_unbounded));
// partial pooling of observations across groups/geographies (exchangeable "random effects")
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }  
// process model (likelihood of the data) + prior degrees of freedom
  if (is_student) {
    nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);
    y ~ student_t(nu[1], f, sigma[has_sigma]);
  } 
  if (is_gaussian)  y ~ normal(f, sigma[has_sigma]);
  if (is_poisson) y_int ~ poisson_log(f);
  if (is_binomial) y_int ~ binomial(trials, f);


