  if (has_re) {
    for (i in 1:n) {
      fitted[i] += alpha_re[id[i]];
   }
  }  
  if (dwx) {
   if (has_me) {
      for (i in 1:dwx) {
     fitted += csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_all[,wx_idx[i]]) * gamma[i];
     }
   } else {
      fitted += WX * gamma;
      }
  } 
  if (dx_all) fitted += x_all * beta;
  if (is_binomial) fitted = inv_logit(fitted);
  if (is_poisson) fitted = exp(fitted);

