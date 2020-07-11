// covariates to model, percentages/proportions on [0,100] and continuous [-inf,inf] variables
  matrix<lower=0,upper=100>[n, dx_me_prop] x_true_prop;
  matrix[n, dx_me_cont] x_true_cont;
// regression parameters
  real intercept;
  vector[dwx] gamma;
  vector[dx_all] beta;
  real<lower=0> nu[is_student]; 
  real<lower=0> sigma[has_sigma];
// for partial pooling across groups/geographies
  vector[n_ids] alpha_re_tilde;
  real<lower=0> alpha_tau[has_re];


