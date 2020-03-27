// if Poisson, offset term is natural log of user-provided offset vector
  vector[n] f;
  vector[dx] beta;
  vector[dev] beta_ev;  
  real error_scale;
  if (has_sigma) {
    error_scale = sigma[1];
      } else {
    error_scale = 1;
      }
  beta_ev = rhs_prior(dev, z, aux1_global, aux2_global, aux1_local, aux2_local, caux, scale_global, slab_scale, error_scale);
  f = offset + intercept + EV * beta_ev; 
  if (dx) {
    f += Q_ast * beta_tilde;
    beta = R_inverse * beta_tilde;
  }
  if (has_re) {
    for (i in 1:n) {
      f[i] += alpha_tau[has_re] * alpha_re_tilde[id[i]];
   }
  }

