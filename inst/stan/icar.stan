functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=1,upper=3> type; // 1=icar, 2=bym, 3=bym2
  int<lower=1> k; // no. of groups
  int group_size[k]; // observational units per group
  int group_idx[n]; // index of observations, ordered by group
  int<lower=1> n_edges; 
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  vector[n_edges] weight;
  int<lower=1, upper=k> comp_id[n]; 
  vector<lower=1>[k] scale_factor;
}

transformed data {
  int has_theta = type > 1;
#include parts/trans_data.stan
}

parameters {
  vector[n] phi_tilde;
  real<lower=0> spatial_scale;
  vector[type > 1 ? n : 0] theta_tilde;
  real<lower=0> theta_scale[type == 2];
  real<lower=0,upper=1> rho[type == 3];  
#include parts/params.stan
}

transformed parameters {
  vector[type > 1 ? n : 0] convolution;
#include parts/trans_params_declaration.stan
  if (type == 1) f += phi_tilde * spatial_scale;
  if (type == 2) {
    convolution = convolve_bym(phi_tilde * spatial_scale, theta_tilde * theta_scale[1], n, k, group_size, group_idx);
    f += convolution;
  }  
  if (type == 3) {
    convolution = convolve_bym2(phi_tilde, theta_tilde, spatial_scale, n, k, group_size, group_idx, rho[1], scale_factor);
    f += convolution;
  }  
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
  if (has_theta) {
   theta_tilde ~ std_normal();
   if (type == 2) theta_scale[1] ~ std_normal();
   if (type == 3) rho[1] ~ beta(1, 1);
  }
  spatial_scale ~ std_normal();
  phi_tilde ~ icar_normal(node1, node2, k, group_size, group_idx, has_theta);  
 }

generated quantities {
  vector[n] phi; 
  vector[type > 1 ? n : 0] theta; 
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
    if (type < 3) phi[i] = phi_tilde[i] * spatial_scale;
    if (type == 2) theta[i] = theta_tilde[i] * theta_scale[1];    
    if (type == 3) {
      phi[i] = spatial_scale * sqrt( rho[1] * inv(scale_factor[comp_id[i]]) ) * phi_tilde[i];
      theta[i] = spatial_scale * sqrt(1 - rho[1]) * theta_tilde[i];
    }   
#include parts/gen_quants_expression_in_loop.stan      
  }
}

