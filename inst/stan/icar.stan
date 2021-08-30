functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=1,upper=3> type; // 1=icar, 2=bym, 3=bym2
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
}

transformed data {
  int has_theta = type > 1;
#include parts/trans_data.stan
}

parameters {
  vector[m] alpha_phi;  
  vector[n] phi_tilde;
  real<lower=0> spatial_scale;
  vector[type > 1 ? n : 0] theta_tilde;
  real<lower=0> theta_scale[type == 2];
  real<lower=0,upper=1> rho[type == 3];  
#include parts/params.stan
}

transformed parameters {
  vector[n] phi;
  vector[type > 1 ? n : 0] theta;
#include parts/trans_params_declaration.stan
  if (type == 1) {
    phi = make_phi(phi_tilde, spatial_scale, 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += phi;
  }
  if (type == 2) {
    theta = theta_tilde * theta_scale[1];
    phi = make_phi(phi_tilde, spatial_scale, 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym(phi, theta, n, k, group_size, group_idx);
  }  
  if (type == 3) {
    theta = spatial_scale * sqrt(1 - rho[1]) * theta_tilde;
    phi = make_phi(phi_tilde, spatial_scale, rho[1], inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym2(phi_tilde, theta_tilde, spatial_scale, n, k, group_size, group_idx, rho[1], inv_sqrt_scale_factor);    
  }  
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
  if (has_theta) {
    target += std_normal_lpdf(theta_tilde);
    if (type == 2) target += std_normal_lpdf(theta_scale[1]);
   // implicit uniform prior on rho:   if (type == 3) rho[1] ~ beta(1, 1);
  }
  target += std_normal_lpdf(spatial_scale);
  phi_tilde ~ icar_normal(spatial_scale, node1, node2, k, group_size, group_idx, has_theta);
  if (m) target += normal_lpdf(alpha_phi | 0, alpha_prior[2]);  
 }

generated quantities {
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan      
  }
}

