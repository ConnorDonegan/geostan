functions {
#include sar-lpdf.stan
}

data {
  // data
  int<lower=1> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n, k] x;


// SAR
  int nW_w;
  vector[nW_w] W_w;
  array[nW_w] int W_v;
  array[n + 1] int W_u;
  vector[n] eigenvalues_w;
}

parameters {
  // SA parameter
  real<lower=1/min(eigenvalues_w), upper=1/max(eigenvalues_w)> rho;     

  // scale parameter
  real<lower=0> sigma;

  // intercept
  real alpha;

  // coefficients
  vector[k] beta;
}

model{
  vector[n] mu = alpha + x * beta;

  // Likelihood: Y ~ Normal(Mu, Sigma)
  target += sar_normal_lpdf(y |
                  mu, sigma, rho,
                  W_w,
                  W_v,
                  W_u,
                  eigenvalues_w,
                  n);

  // prior for scale parameter
  target += student_t_lpdf(sigma | 10, 0, 5);

  // prior for beta
  target += normal_lpdf(beta | 0, 5);

  // prior for intercept
  target += normal_lpdf(alpha | 0, 5);
}
