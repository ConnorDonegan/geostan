// likelihood function
  int<lower=0,upper=6> family;

// number of observations
  int<lower=0> n; 

  int<lower=0> center_x;

// missing y
  int<lower=0> n_mis;
  int<lower=0> n_obs;
  array[n_mis] int y_mis_idx;
  array[n_obs] int y_obs_idx;

// censored counts
  int censor_point;

// outcome 
  vector[n] y;
  array[n] int<lower=0> y_int;
  array[n] int<lower=0> trials;
  int<lower=0,upper=1> prior_only;

// offset
  vector[n] input_offset; 

  // tells which of X are in SLX
  int<lower=0> dwx;
  array[dwx ? dwx : 1] int wx_idx; 

  // connectivity matrix: (row-standardized) for spatial lag of X or SAR
  int<lower=0> nW_w;
  vector[nW_w] W_w;
  array[nW_w] int W_v;
  array[nW_w > 1 ? n + 1 : 1] int W_u;
  
  // covariates and observational error information
    // lower, upper bounds for bounded data models
  array[2] real bounds;
    // no. columns
  int<lower=0> dx_obs;
  int<lower=0> dx_me;
  array[dx_me] int<lower=0,upper=1> use_logit;
    // indices matching columns of observed and ME data matrices to columns of raw data matrix x (and parameter x_all)
  array[dx_obs ? dx_obs : 1] int<lower=0> x_obs_idx; 
  array[dx_me ? dx_me : 1] int<lower=0> x_me_idx;  
    // covariates observed with practical certainty 
  matrix[n, dx_obs ? dx_obs : 0] x_obs;
    // covariates observed with uncertainty, and standard errors
  array[dx_me] vector[n] x_me;
  array[dx_me] vector<lower=0>[n] sigma_me;
  //  priors for x_true
  vector[dx_me] prior_nux_true_alpha;
  vector[dx_me] prior_nux_true_beta;
  vector[dx_me] prior_mux_true_location;
  vector[dx_me] prior_mux_true_scale;
  vector[dx_me] prior_sigmax_true_df;
  vector[dx_me] prior_sigmax_true_location;
  vector[dx_me] prior_sigmax_true_scale;
  vector[2] prior_rhox_true;

  // data for auto-normal [ME] models
  int<lower=0,upper=1> spatial_me;
  int<lower=0,upper=1> WCAR;
  int nA_w;
  vector[nA_w] A_w;
  array[nA_w] int A_v;
  array[n + 1] int A_u;
  vector[n] Delta_inv;
  real log_det_Delta_inv;
  vector[n] lambda;

  // non-spatial partial pooling 
  int<lower=0,upper=1> has_re; // has varying intercept?
  int<lower=0> n_ids; // number of units
  array[n] int<lower=0,upper=n_ids> id; // identifier for the observational units associated with the varying intercepts

  // priors
  vector[2] prior_alpha; // prior on the intercept
  int<lower=0> dbeta_prior;
  vector[dbeta_prior] prior_beta_location; // coefficient priors, with any SLX terms listed first; 
  vector<lower=0>[dbeta_prior] prior_beta_scale;
  //row_vector[dbeta_prior] beta_prior[2]; 
  vector[3] prior_alpha_tau; // prior on standard deviation of varying intercepts
  vector[2] prior_t_nu;
  vector[3] prior_sigma;

  // ICAR 
  int<lower=0,upper=3> type; // 0=glm, 1=icar, 2=bym, 3=bym2
  int<lower=1> k; // no. of groups
  array[k] int group_size; // observational units per group
  array[n] int group_idx; // index of observations, ordered by group
  int<lower=0> m; // no of components requiring additional intercepts
  matrix[n, m] A; // dummy variables for any extra graph component intercepts  
  int<lower=1> n_edges; 
  array[n_edges] int<lower=1, upper=n> node1;
  array[n_edges] int<lower=1, upper=n> node2;
  vector[n_edges] weight;
  array[n] int<lower=1, upper=k> comp_id; 
  vector[k] inv_sqrt_scale_factor; // can be a vector of ones, as a placeholder

  // ESF 
  int<lower=0> dev; // number of eigenvectors : now included in parts/data.stan
  matrix[n, dev] EV; // the eigenvectors : now included in parts/data.stan
  real<lower=0> global_scale;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df;

  // CAR
  array[2] real car_rho_lims;
  int<lower=0,upper=1> car;

  // CAR or SAR: zero-mean parameterized
  int<lower=0, upper=1> ZMP;

  // SAR: use W matrix (as above: W_w, W_u, W_v)
//  int nImW_w;
//  int nW;
//  vector[nImW_w] ImW_w;
//  array[nImW_w] int ImW_v;
//  array[n + 1] int ImW_u;
//  array[nW] int Widx;
  vector[n] eigenvalues_w;
  array[2] real sar_rho_lims;
  int<lower=0,upper=1> sar;

