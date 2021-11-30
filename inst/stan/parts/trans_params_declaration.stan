  // ESF
  vector[dev] beta_ev;
  vector[dev ? n : 0] esf;  
  real error_scale[dev ? 1 : 0];  
  // ICAR
  vector[type ? n : 0] phi;
  vector[type > 1 ? n : 0] theta;
  // CAR
  vector[car ? n : 0] log_lambda_mu;  
  // GLM
  matrix[n, dx_all] x_all;
  vector[n] fitted;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me) for (j in 1:dx_me) x_all[ ,x_me_idx[j]] = x_true[j];
  if (center_x) for (j in 1:dx_all) x_all[,j] = x_all[,j] - mean(x_all[,j]);
  if (!car) fitted = offset + intercept;  



