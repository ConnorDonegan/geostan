functions {
#include wcar-lpdf.stan
}

data {
  // data
  int<lower=1> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n, k] x;

    // CAR parts
  int nA_w;
  vector[nA_w] A_w;
  array[nA_w] int A_v;
  array[n + 1] int A_u;
  vector[n] Delta_inv;
  real log_det_Delta_inv;
  vector[n] lambda;
}

parameters {
  // spatial autocorrelation (SA) parameter
  real<lower=1/min(lambda), upper=1/max(lambda)> rho;
  
  // scale parameter
  real<lower=0> tau;
  
  // intercept
  real alpha;
  
  // coefficients
  vector[k] beta;  
}

model {
  vector[n] mu = alpha + x * beta;

  // Likelihood: y ~ Normal(Mu, Sigma)
  target += wcar_normal_lpdf(y |
                 mu, tau, rho, // mean, scale, SA
                 A_w, A_v, A_u, // stuff from prep_car_data
                 Delta_inv, 
                 log_det_Delta_inv,
                 lambda, 
                 n);    
  
  // prior for scale parameter
  target += student_t_lpdf(tau | 10, 0, 5);

  // prior for beta
  target += normal_lpdf(beta | 0, 5);

  // prior for intercept
  target += normal_lpdf(alpha | 0, 5);
}
