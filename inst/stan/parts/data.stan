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
  real bounds[2];
    // no. columns
  int<lower=0> dx_obs;
  int<lower=0> dx_me;
    // indices matching columns of observed and ME data matrices to columns of raw data matrix x (and parameter x_all)
  int<lower=0> x_obs_idx[dx_obs ? dx_obs : 1]; 
  int<lower=0> x_me_idx[dx_me ? dx_me : 1];  
    // covariates observed with practical certainty 
  matrix[n, dx_obs ? dx_obs : 0] x_obs;
    // covariates observed with uncertainty, and standard errors
  vector[n] x_me[dx_me];
  vector<lower=0>[n] sigma_me[dx_me];
  //  priors for x_true
  vector[dx_me] prior_nux_true_alpha;
  vector[dx_me] prior_nux_true_beta;
  vector[dx_me] prior_mux_true_location;
  vector[dx_me] prior_mux_true_scale;
  vector[dx_me] prior_sigmax_true_df;
  vector[dx_me] prior_sigmax_true_location;
  vector[dx_me] prior_sigmax_true_scale;
  vector[2] prior_rhox_true;
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
  vector[2] prior_alpha; // prior on the intercept
  int<lower=0> dbeta_prior;
  vector[dbeta_prior] prior_beta_location; // coefficient priors, with any SLX terms listed first; 
  vector<lower=0>[dbeta_prior] prior_beta_scale;
//row_vector[dbeta_prior] beta_prior[2]; 
  vector[3] prior_alpha_tau; // prior on standard deviation of varying intercepts
  vector[2] prior_t_nu;
  vector[3] prior_sigma;
  // outcome 
  vector[n] y;
  int<lower=0> y_int[n];
  int<lower=0> trials[n];
  int<lower=0,upper=1> prior_only;

  // ICAR 
  int<lower=0,upper=3> type; // 0=glm, 1=icar, 2=bym, 3=bym2
  int<lower=1> k; // no. of groups
  int group_size[k]; // observational units per group
  int group_idx[n]; // index of observations, ordered by group
  int<lower=0> m; // no of components requiring additional intercepts
  matrix[n, m] A; // dummy variables for any extra graph component intercepts  
  int<lower=1> n_edges; 
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  vector[n_edges] weight;
  int<lower=1, upper=k> comp_id[n]; 
  vector[k] inv_sqrt_scale_factor; // can be a vector of ones, as a placeholder

  // ESF 
  int<lower=0> dev; // number of eigenvectors : now included in parts/data.stan
  matrix[n, dev] EV; // the eigenvectors : now included in parts/data.stan
  real<lower=0> global_scale;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df;

  // CAR
  real car_rho_lims[2];
  int<lower=0,upper=1> car;
