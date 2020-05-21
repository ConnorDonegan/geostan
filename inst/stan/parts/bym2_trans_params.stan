// if Poisson, offset term is natural log of user-provided offset vector
  vector[dx] beta;
  real<lower=0, upper=1> rho;
  vector[n] convolved_re;
  vector[n] f;
  rho = inv_logit(logit_rho);
  convolved_re = sigma_re * (sqrt(rho / scaling_factor) * v + sqrt(1 - rho) * u);
  f = offset + intercept + convolved_re; 
  if (dx) {
    f += Q_ast * beta_tilde;
    beta = R_inverse * beta_tilde;
  }
  if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
   }
  }

