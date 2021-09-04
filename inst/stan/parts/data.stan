// likelihood function
  int<lower=0,upper=5> family;
// number of observations
  int<lower=0> n; 
  
// offset
  vector[n] offset; 
  
// connectivity matrix: row-standardized for spatial lag of X 
  int<lower=0> dwx;
  int wx_idx[dwx ? dwx : 1];
  int<lower=0> dw_nonzero;
  vector[dw_nonzero] W_w;
  int W_v[dw_nonzero];
  int W_u[dwx ? n + 1 : 1];
  
// covariates and observational error information
    // lower, upper bounds for bounded data models
  vector[2] bounds;
    // no. columns
  int<lower=0> dx_obs;
  int<lower=0> dx_me_bounded;
  int<lower=0> dx_me_unbounded;
    // indices matching columns of observed and ME data matrices to columns of raw data matrix x (and parameter x_all)
  int<lower=0> x_obs_idx[dx_obs ? dx_obs : 1];
  int<lower=0> x_me_unbounded_idx[dx_me_unbounded ? dx_me_unbounded : 1];
  int<lower=0> x_me_bounded_idx[dx_me_bounded ? dx_me_bounded : 1];
    // covariates observed with practical certainty 
  matrix[n, dx_obs ? dx_obs : 0] x_obs;
    // covariates observed with uncertainty, and standard errors: unbounded
  vector[n] x_me_unbounded[dx_me_unbounded];
  vector<lower=0>[n] sigma_me_unbounded[dx_me_unbounded];
   // covariates observed with uncertainty, and standard errors: bounded
  vector<lower=bounds[1],upper=bounds[2]>[n] x_me_bounded[dx_me_bounded];
  vector<lower=0>[n] sigma_me_bounded[dx_me_bounded];
  //  priors for x_true 
  vector[dx_me_unbounded] prior_mux_true_unbounded_location;
  vector[dx_me_unbounded] prior_mux_true_unbounded_scale;
  vector[dx_me_unbounded] prior_sigmax_true_unbounded_scale;
  vector[dx_me_bounded] prior_mux_true_bounded_location;
  vector[dx_me_bounded] prior_mux_true_bounded_scale;
  vector[dx_me_bounded] prior_sigmax_true_bounded_scale;

  // data for auto-Guassian [ME] models
  int<lower=0,upper=1> spatial_me;
  int<lower=0,upper=1> WCAR;
  int nAx_w;
  int nC;
  vector[nAx_w] Ax_w;
  int Ax_v[nAx_w];
  int Ax_u[n + 1];
  int Cidx[nC];
  vector[n] Delta_inv;
  real log_det_Delta_inv;
  vector[n] lambda;

// non-spatial partial pooling 
  int<lower=0,upper=1> has_re; // has varying intercept?
  int<lower=0> n_ids; // number of units
  int<lower=0,upper=n_ids> id[n]; // identifier for the observational units associated with the varying intercepts
// priors
  vector[2] alpha_prior; // prior on the intercept
  int<lower=0> dbeta_prior;  
  row_vector[dbeta_prior] beta_prior[2]; // coefficient priors, with any SLX terms listed first; 
  vector[3] alpha_tau_prior; // prior on standard deviation of varying intercepts
  vector[2] t_nu_prior;
  vector[3] sigma_prior;
  // outcome 
  vector[n] y;
  int<lower=0> y_int[n];
  int<lower=0> trials[n];
  int<lower=0,upper=1> prior_only;

