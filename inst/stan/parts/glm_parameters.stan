  real intercept;
  vector[dx] beta_tilde;
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];
  real<lower=0> nu[is_student];
  real<lower=0> sigma[has_sigma];

