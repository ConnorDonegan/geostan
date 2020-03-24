// esf data
  int<lower=0> dev; // number of eigenvectors 
  matrix[n, dev] EV; // the eigenvectors
  real<lower=0> scale_global;  // horseshoe parameters
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
