  intercept ~ student_t(alpha_prior[1], alpha_prior[2], alpha_prior[3]);
  if (dx) beta ~ student_t(beta_prior[1], beta_prior[2], beta_prior[3]);
  if (has_re) {
    alpha_tau[has_re] ~ student_t(alpha_tau_prior[1], alpha_tau_prior[2], alpha_tau_prior[3]);
    alpha_re_tilde ~ std_normal();    
  }
