// any other declarations in this section must occur *before* this file
  matrix[n, dx_all] x_all;
  vector[n] f;
  if (dx_obs) x_all[,x_obs_idx] = x_obs;
  if (dx_me_unbounded) for (j in 1:dx_me_unbounded) x_all[ ,x_me_unbounded_idx[j]] = x_true_unbounded[j];
  if (dx_me_bounded) for (j in 1:dx_me_bounded) x_all[,x_me_bounded_idx[j]] = x_true_bounded[j];
  if (model_offset) {
     f = offset_est + intercept;
     } else {
     f = offset_obs + intercept;  // default offset_obs=0
   }

