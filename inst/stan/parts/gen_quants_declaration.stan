 // any other declarations in this block must be made *before* this included file
  vector[n] log_lik;
  vector[n] yrep;
  vector[n] residual;
  vector[n] fitted;
  vector[n_ids] alpha_re;
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }




