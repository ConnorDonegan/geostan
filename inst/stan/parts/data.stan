// likelihood function
  int<lower=0,upper=4> family;
// number of observations
  int<lower=0> n; 
  
// offset
  vector[n] offset; 
  
// connectivity matrix
  int<lower=0> dwx;
  int wx_idx[dwx ? dwx : 1];
  int<lower=0> dw_nonzero;
  matrix[dwx ? n : 1, dwx ? n : 1] W;
  
// covariates and measurement error information
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
    // `observed' data and ME data (unbounded and bounded)
  matrix[n, dx_obs ? dx_obs : 0] x_obs;
  vector[n] x_me_unbounded[dx_me_unbounded];
  vector<lower=bounds[1],upper=bounds[2]>[n] x_me_bounded[dx_me_bounded];
    // standard errors of ME data
  vector<lower=0>[n] sigma_me_unbounded[dx_me_unbounded];
  vector<lower=0>[n] sigma_me_bounded[dx_me_bounded];
  // Spatial observational error model (ME)
  int<lower=1> dev; // lower=1 to allow for an empty EV matrix (is column vector of zeros when not needed)
  matrix[n,dev] EV; // eigenvectors for spatial filter
  int<lower=0,upper=1> spatial_me;
    // priors for unbounded variables
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] slab_df_me_unbounded;  
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] slab_scale_me_unbounded;
  vector<lower=0>[spatial_me ? dx_me_unbounded : 0] scale_global_me_unbounded;
   // priors for bounded variables
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] slab_df_me_bounded;    
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] slab_scale_me_bounded;  
  vector<lower=0>[spatial_me ? dx_me_bounded : 0] scale_global_me_bounded;
  
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

