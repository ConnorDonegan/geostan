 // any other declarations in this section must occur *before* this file
  matrix[n, dx_all] x_all;
  vector[n] f;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me_bounded) x_all[,x_me_bounded_idx] = x_true_bounded;
  if (dx_me_unbounded) x_all[,x_me_unbounded_idx] = x_true_unbounded;
  if (!model_offset && !is_poisson) f = offset_obs;     // default offset_obs=0
  if (model_offset && !is_poisson) f = offset_est;
  if (is_poisson) {
     if (model_offset) f = log(offset_est);
     if (!model_offset && has_offset) f = log(offset_obs);
  }
  f += intercept; 
  

