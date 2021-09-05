functions {
#include parts/priors.stan
}

data {
  //GLM 
#include parts/data.stan
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
  real<lower=0> scale_global;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
}

transformed data {
  int has_theta = type > 1;
#include parts/trans_data.stan
}

parameters {
  //ICAR
  vector[m] alpha_phi;  
  vector[type ? n : 0] phi_tilde;
  real<lower=0> spatial_scale[type ? 1 : 0];
  vector[type > 1 ? n : 0] theta_tilde;
  real<lower=0> theta_scale[type == 2];
  real<lower=0,upper=1> rho[type == 3];
  // ESF
  real<lower=0> aux1_global[dev ? 1 : 0];
  real<lower=0> aux2_global[dev ? 1 : 0];
  vector<lower=0>[dev] aux1_local;
  vector<lower=0>[dev] aux2_local;
  real<lower=0> caux[dev ? 1 : 0];
  vector[dev] z;
  // GLM
#include parts/params.stan
}

transformed parameters {
  // ESF
  vector[dev] beta_ev;  
  real error_scale[dev ? 1 : 0];
  // ICAR
  vector[type ? n : 0] phi;
  vector[type > 1 ? n : 0] theta;
  // GLM
#include parts/trans_params_declaration.stan
  // ICAR
  if (type == 1) {
    phi = make_phi(phi_tilde, spatial_scale[1], 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += phi;
  }
  if (type == 2) {
    theta = theta_tilde * theta_scale[1];
    phi = make_phi(phi_tilde, spatial_scale[1], 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym(phi, theta, n, k, group_size, group_idx);
  }  
  if (type == 3) {
    theta = spatial_scale[1] * sqrt(1 - rho[1]) * theta_tilde;
    phi = make_phi(phi_tilde, spatial_scale[1], rho[1], inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym2(phi_tilde, theta_tilde, spatial_scale[1], n, k, group_size, group_idx, rho[1], inv_sqrt_scale_factor);    
  }
  // ESF
  if (dev) {
    if (has_sigma) {
      error_scale[1] = sigma[1];
    } else {
      error_scale[1] = 1;
    }
    beta_ev = rhs_prior(dev, z, aux1_global[1], aux2_global[1], aux1_local, aux2_local, caux[1], scale_global, slab_scale, error_scale[1]);
  fitted += EV * beta_ev;
  }
#include parts/trans_params_expression.stan
}

model {
  //GLM
#include parts/model.stan
  // ICAR
  if (type) {
    if (has_theta) {
      target += std_normal_lpdf(theta_tilde);
      if (type == 2) target += std_normal_lpdf(theta_scale[1]);
      // implicit uniform prior on rho:   if (type == 3) rho[1] ~ beta(1, 1);
    }
    target += std_normal_lpdf(spatial_scale[1]);
    phi_tilde ~ icar_normal(spatial_scale[1], node1, node2, k, group_size, group_idx, has_theta);
    if (m) target += normal_lpdf(alpha_phi | 0, alpha_prior[2]);
  }
  // ESF
  if (dev) {
    target += std_normal_lpdf(z);
    target += std_normal_lpdf(aux1_local);
    target += inv_gamma_lpdf(aux2_local | 0.5, 0.5); // .5 * nu_local, .5 * nu_local, nu_local = 1
    target += std_normal_lpdf(aux1_global[1]);
    target += inv_gamma_lpdf(aux2_global[1] | 0.5, 0.5); // .5 * nu_local, .5 * nu_global, both = 1
    target += inv_gamma_lpdf(caux[1] | 0.5*slab_df, 0.5*slab_df);
  }
}

generated quantities {
  vector[dev ? n : 0] esf;
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
    if (dev) esf[i] = EV[i] * beta_ev;
#include parts/gen_quants_expression_in_loop.stan      
  }
}

