// any other declarations in this section must occur *before* this file
  matrix[n, dx_all] x_all;
  vector[n] fitted;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me) for (j in 1:dx_me) x_all[ ,x_me_idx[j]] = x_true[j];
  fitted = offset + intercept;

