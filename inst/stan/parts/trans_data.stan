  // ICAR
  int has_theta = type > 1;

  // outcome model information
  int<lower=0,upper=1> is_gaussian;
  int<lower=0,upper=1> is_student;
  int<lower=0,upper=1> is_poisson;
  int<lower=0,upper=1> is_binomial;
  int<lower=0,upper=1> is_auto_gaussian;
  int<lower=0,upper=1> has_sigma;
  int<lower=0,upper=1> has_offset;
  int<lower=0> dx_all;
  int<lower=0> has_me;
  is_gaussian = family == 1;
  is_student =  family == 2;
  is_poisson =  family == 3;
  is_binomial = family == 4;
  is_auto_gaussian = family == 5 || family == 6;  
  has_sigma  =  family < 3;
  has_offset = sum(offset) != 0;
  dx_all = dx_obs + dx_me;
  has_me = dx_me > 0;

