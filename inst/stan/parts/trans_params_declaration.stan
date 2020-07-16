 // any other declarations in this section must occur *before* this file
  matrix[n, dx_all] x_all;
  vector[n] f;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me_prop) x_all[,x_me_prop_idx] = x_true_prop;
  if (dx_me_cont) x_all[,x_me_cont_idx] = x_true_cont;
  if (!model_offset && !is_poisson) f = offset_obs;     // default offset_obs=0
  if (model_offset && !is_poisson) f = offset_est;
  if (is_poisson) {
     if (model_offset) f = log(offset_est);
     if (!model_offset && has_offset) f = log(offset_obs);
  }
  f += intercept; 
  

