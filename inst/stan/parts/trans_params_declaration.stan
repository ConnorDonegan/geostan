  // ESF
  vector[dev] beta_ev;
  vector[dev ? n : 0] esf;  
  array[dev ? 1 : 0] real error_scale;

  // ICAR
  vector[type ? n : 0] phi;
  vector[type > 1 ? n : 0] theta;

  // Hierarchical CAR model with modeled mean (i.e., not constrained to have mean of zero)
  vector[(car||sar) && is_auto_gaussian == 0 && ZMP == 0 ? n : 0] log_lambda_mu;

  // GLM
  matrix[n, use_qr ? 0 : dx_all] x_all;
  vector[n] fitted;
  vector[dx_all] beta;
  vector[dwx] gamma;
  if (use_qr == 0) {
    if (dx_obs) x_all[,x_obs_idx] = x_obs;
    if (dx_me) for (j in 1:dx_me) x_all[ ,x_me_idx[j]] = x_true[j];
    if (center_x) for (j in 1:dx_all) x_all[,j] = x_all[,j] - mean(x_all[,j]);
  }


