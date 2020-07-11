 // any other declarations in this section must occur *before* this file
  matrix[n, dx_all] x_all;
  vector[n] offset;
  vector[n] f;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me_prop) x_all[,x_me_prop_idx] = x_true_prop;
  if (dx_me_cont) x_all[,x_me_cont_idx] = x_true_cont;
  if (!model_offset) offset = offset_obs;
  if (is_poisson && has_offset) offset = log(offset);
  f = offset + intercept; // default offset = 0



