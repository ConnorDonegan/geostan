  // ICAR
  int has_theta = type > 1;

  // QR transformation 
  int<lower=0,upper=1> use_qr;
  int<lower=0> d_qr = dx_obs + dwx;
  matrix[n, d_qr] Q_ast;
  matrix[d_qr, d_qr] R_ast;
  matrix[d_qr, d_qr] R_ast_inverse;
  matrix[n, d_qr] x_tmp; // [X WX]

  // outcome model information
  int<lower=0,upper=1> is_gaussian;
  int<lower=0,upper=1> is_student;
  int<lower=0,upper=1> is_poisson;
  int<lower=0,upper=1> is_binomial;
  int<lower=0,upper=1> is_auto_gaussian;
  int<lower=0,upper=1> has_sigma;
  int<lower=0,upper=1> has_offset;
  int<lower=0> dx_all;
  int<lower=0,upper=1> has_me;
  is_gaussian = family == 1;
  is_student =  family == 2;
  is_poisson =  family == 3;
  is_binomial = family == 4;
  is_auto_gaussian = family == 5 || family == 6;  
  has_sigma  =  family < 3;
  has_offset = sum(input_offset) != 0;
  dx_all = dx_obs + dx_me;
  has_me = dx_me > 0;
  use_qr = dx_obs > 0 && has_me == 0;

  if (use_qr) {
   	
    if (center_x) { // center covariates
      for (j in 1:dx_obs) x_tmp[,j] = x_obs[,j] - mean(x_obs[,j]);
    } else {
      x_tmp[ , 1:dx_obs ] = x_obs;
    }

    if (dwx) { // add spatially-lagged covariates
      for (j in 1:dwx) x_tmp[ , dx_obs + j ] = csr_matrix_times_vector(n, n, W_w, W_v, W_u, x_obs[,wx_idx[j]]);       
    }
    
    // create, thin, and scale the QR decomposition
    Q_ast = qr_thin_Q(x_tmp) * sqrt(n - 1);
    R_ast = qr_thin_R(x_tmp) / sqrt(n - 1);
    R_ast_inverse = inverse(R_ast);
  }
