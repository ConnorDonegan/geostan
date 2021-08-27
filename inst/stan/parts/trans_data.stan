// outcome model information
  int<lower=0,upper=1> is_gaussian;
  int<lower=0,upper=1> is_student;
  int<lower=0,upper=1> is_poisson;
  int<lower=0,upper=1> is_binomial;
  int<lower=0,upper=1> is_auto_gaussian;
  int<lower=0,upper=1> has_sigma;
  int<lower=0,upper=1> has_offset;
  int<lower=0> dx_all;
// get summary data on ME covariates for very weakly informative priors
  int<lower=0,upper=1> has_me;
  vector[dx_me_unbounded] prior_scale_x_true_unbounded;
  vector[dx_me_unbounded] prior_mean_x_true_unbounded;
  vector[dx_me_bounded] prior_scale_x_true_bounded;
  vector[dx_me_bounded] prior_mean_x_true_bounded;
// for WX, sparse matrix representation
  vector[dw_nonzero] w;
  int v[dw_nonzero];
  int u[n + 1];
  matrix[n, dwx] WX;
// for auto-Gaussian models
  vector[n] mean_zero = rep_vector(0, n);

  is_gaussian = family == 1;
  is_student =  family == 2;
  is_poisson =  family == 3;
  is_binomial = family == 4;
  is_auto_gaussian = family == 5;  
  has_sigma  =  family < 3;
  has_offset = sum(offset) != 0;
  dx_all = dx_obs + dx_me_bounded + dx_me_unbounded;
  has_me = (dx_all > dx_obs);
  if (dwx) {
    w = csr_extract_w(W);
    v = csr_extract_v(W);
    u = csr_extract_u(W);
  }
  if ((!has_me) && dwx) {
     WX = W * x_obs[,wx_idx];
  }
  if (dx_me_unbounded) {
    for (j in 1:dx_me_unbounded) {
       prior_scale_x_true_unbounded[j] = sd(x_me_unbounded[j]);
       prior_mean_x_true_unbounded[j] = mean(x_me_unbounded[j]);
    }
  }
  if (dx_me_bounded) {
    for (j in 1:dx_me_bounded) {
       prior_scale_x_true_bounded[j] = sd(x_me_bounded[j]);
       prior_mean_x_true_bounded[j] = mean(x_me_bounded[j]);
    }
  }
  
