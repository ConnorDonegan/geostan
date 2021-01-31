/**
  * Gather terms for the regularized horseshoe model (Piironen and Vehtari)
  *
  * @return A vector of coefficientes
  **/
  vector rhs_prior(int dev,
		   vector z,
		   real aux1_global,
		   real aux2_global,
		   vector aux1_local,
		   vector aux2_local,
		   real caux,
		   real scale_global,
		   real slab_scale,
		   real error_scale) {
  real tau = aux1_global * sqrt(aux2_global) * scale_global * error_scale;
  real c = slab_scale * sqrt(caux);
  vector[dev] lambda = aux1_local .* sqrt(aux2_local);
  vector[dev] lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + square(tau) * square(lambda)) );
  return z .* lambda_tilde * tau;
}

/**
  * Return the log probability of the conditional autoregressive (CAR) prior,
  *  excluding additive constants.
  *
  * @param y Vector containing the parameters with a CAR prior
  * @param mu Mean vector.
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param w sparse representation of W' (transpose!): non-zero values
  * @param v column indices for values in w
  * @param u row starting indices for values in w followed by size of w
  * @param D_diag Diagonal of D matrix; e.g., number of neighbors for each location
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  * @param n Length of y 
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real car_normal_lpdf(vector y, vector mu,
                       real tau, real alpha,
                       vector w, int[] v, int[] u, 
                       vector D_diag, vector lambda,
                       int n) {
  vector[n] yc = y - mu; 
  real ytDy; // yc transpose * D * yc
  real ytWy; // yc transpose * W
  vector[n] ldet_prec;    
  ytDy = (yc .* D_diag)' * yc;
  ytWy = csr_matrix_times_vector(n, n, w, v, u, yc)' * yc;    
  for (i in 1:n) ldet_prec[i] = log1m(alpha * lambda[i]);
  return 0.5 * (n * log(tau)
		+ sum(ldet_prec)
		- tau * (ytDy - alpha * ytWy));
}

/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
 * excluding additive constants. 
 *
 * @param phi Vector of parameters for spatial smoothing (on unit scale, approximately)
 * @param node1 
 * @param node2
 * @param k number of groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param has_theta If the model contains an independent partial pooling term, phi for singletons can be zeroed out; otherwise, those observations require a standard normal prior. Both BYM and BYM2 have theta.
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, 
		      int[] node1, int[] node2, 
		      int k, int[] group_size, int[] group_idx,
		      int has_theta) {
  real lp;
  int pos=1;
  lp = -0.5 * dot_self(phi[node1] - phi[node2]);
  if (has_theta) {
    for (j in 1:k) {
      /* sum to zero constraint for each connected group; singletons zero out */
      lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * group_size[j]);
      pos += group_size[j];
    }
  } else {
    /* has no theta */
    for (j in 1:k) {
      if (group_size[j] > 1) {
    /* same as above for non-singletons: sum to zero constraint */
    lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * group_size[j]);
      } else {
    /* its a singleton: independent std normal prior */
    lp += std_normal_lpdf(phi[segment(group_idx, pos, group_size[j])]);
      }      
      pos += group_size[j];
    }
  }
  return lp;
}


/**
 * Combine local and global partial-pooling components into the convolved BYM term.
 *
 * @param phi local (spatially autocorrelated) component
 * @param theta global component
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 *
 * @return BYM convolution vector
 */
vector convolve_bym(vector phi, vector theta,
		    int n, int k,
		    int[] group_size, int[] group_idx
		      ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
     if (group_size[j] == 1) {
        convolution[ segment(group_idx, pos, group_size[j]) ] = theta[ segment(group_idx, pos, group_size[j]) ];
    } else {
    convolution[ segment(group_idx, pos, group_size[j]) ] =
     phi[ segment(group_idx, pos, group_size[j]) ] + theta[ segment(group_idx, pos, group_size[j]) ];
  }
      pos += group_size[j];
  }
  return convolution;
}


/**
 * Combine local and global partial-pooling components into the convolved BYM2 term.
 *
 * @param phi_tilde local (spatially autocorrelated) component
 * @param theta_tilde global component
 * @param spatial_scale scale parameter for the convolution term
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param rho proportion of convolution that is spatially autocorrelated
 * @param scale_factor The scaling factor for the BYM2 model. See scale_c R function, using R-INLA.
 *
 * @return BYM2 convolution vector
 */
vector convolve_bym2(vector phi_tilde, vector theta_tilde,
		     real spatial_scale,
		     int n, int k,
		     int[] group_size, int[] group_idx,
		     real rho, vector scale_factor
		      ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
    if (group_size[j] == 1) {
        convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * theta_tilde[ segment(group_idx, pos, group_size[j]) ];
    } else {
    convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * (
     sqrt(rho / scale_factor[j]) * phi_tilde[ segment(group_idx, pos, group_size[j]) ] +
     sqrt(1 - rho) * theta_tilde[ segment(group_idx, pos, group_size[j]) ]
      );
  }
  pos += group_size[j];
  }
  return convolution;
}

