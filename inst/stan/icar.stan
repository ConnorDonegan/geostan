functions {
#include parts/priors.stan
}

data {
#include parts/data.stan
  int<lower=1,upper=3> type; // 1=iar, 2=bym, 3=bym2
  int<lower=1> k; // no. of groups
  int group_size[k]; // observational units per group
  int group_idx[n]; // index of observations, ordered by group
  int<lower=0> n_edges; 
  int<lower=1, upper=n> node1[n_edges];
  int<lower=1, upper=n> node2[n_edges];
  vector<lower=0>[n] scale_factor;
  real<lower=0> phi_scale_prior; 
}

transformed data {
#include parts/trans_data.stan
}

parameters {
//  vector[n] phi_tilde;
//  real<lower=0> phi_scale;
  vector[n] phi;
  real<lower=0> spatial_scale;
  vector[type > 1 ? n : 0] theta;
  real<lower=0> theta_scale[type == 2 ? 1 : 0];
  real logit_rho[type == 3 ? 1 : 0];
#include parts/params.stan
}

transformed parameters {
  real<lower=0, upper=1> rho[type == 3 ? 1 : 0];
#include parts/trans_params_declaration.stan
  if (type == 1) f += phi * spatial_scale;
  if (type == 2) f += phi * spatial_scale + theta * theta_scale[1];
  if (type == 3) {
    rho[1] = inv_logit(logit_rho[1]);
    f += ( sqrt(rho[1] * inv(scale_factor)) .* phi + sqrt(1 - rho[1]) * theta ) * spatial_scale;
  }
//  f += phi_tilde * phi_scale;
#include parts/trans_params_expression.stan
}

model {
  int pos;
#include parts/model.stan
  pos = 1;
  if (type == 1) {
    for (j in 1:k) {
      if (group_size[j] > 1) {
        sum(phi[segment(group_idx, pos, group_size[j])]) ~ normal(0, 0.001 * group_size[j]);
	} else {
	phi[segment(group_idx, pos, group_size[j])] ~ std_normal();
        }
      pos = pos + group_size[j];
    }
  }
  if (type > 1) {
    for (j in 1:k) {
      sum(phi[segment(group_idx, pos, group_size[j])]) ~ normal(0, 0.001 * group_size[j]);
      pos = pos + group_size[j];
    }
   theta ~ std_normal();
   if (type == 2) theta_scale[1] ~ std_normal();
   if (type == 3) logit_rho[1] ~ std_normal();
   }
   target += -0.5 * dot_self(phi[node1] - phi[node2]);
   spatial_scale ~ std_normal();
// IAR model
//  phi_tilde ~ icar_normal(n, node1, node2);
//  phi_scale ~ normal(0, phi_scale_prior);
 }

generated quantities {
//  vector[n] phi;
  vector[n] ssre; // spatially structured partial-pooling term
  vector[type > 1 ? n : 0] sure; // unstructure partial-pooling term
  vector[type > 1 ? n : 0] convolved_re; // ssre + sure
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
 //     phi[i] = phi_tilde[i] * phi_scale;
    if (type == 1) ssre[i] = phi[i] * spatial_scale;
    if (type == 2) {
      ssre[i] = phi[i] * spatial_scale;
      sure[i] = theta[i] * theta_scale[1];
      convolved_re[i] = ssre[i] + sure[i];
    }
    if (type == 3) {
      ssre[i] = spatial_scale * sqrt(rho[1] * inv(scale_factor[i])) * phi[i];
      sure[i] = spatial_scale * sqrt(1 - rho[1]) * theta[i];
      convolved_re[i] = ssre[i] + sure[i];      
    }   
#include parts/gen_quants_expression_in_loop.stan      
  }
}

