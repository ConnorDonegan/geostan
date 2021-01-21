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
 * Log probability of a conditional autoregressive (CAR) model,
 * dropping additive constants.
 *
 *           y ~ N(mu, tau * (D - alph * W)^(-1))
 *
 * @param y Vector containing the parameters with a CAR prior
 * @param mu Mean vector.
 * @param tau Precision parameter for the CAR prior (real)
 * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
 * @param w sparse representation of W' (transpose!): contains all non-zero values of W.
 * @param v column indices for values in w
 * @param u row starting indices for values in w followed by size of w
 * @param D_diag Diagonal of D matrix; e.g., number of neighbors for each location
 * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
 * @param n Length of y 
 *
 * @return Log probability density of CAR model up to additive constant
 */
real car_normal_lpdf(vector y, vector mu,
		     real tau, real alpha,
		     vector w, int[] v, int[] u, 
		     vector D_diag, vector lambda,
		     int n) {
  vector[n] yc = y - mu; 
  row_vector[n] yct_D; // yc transpose * D
  row_vector[n] yct_W; // yc transpose * W
  vector[n] ldet_terms;    
  yct_D = (yc .* D_diag)';
  yct_W = csr_matrix_times_vector(n, n, w, v, u, yc)';    
  for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
  return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
		- tau * (yct_D * yc - alpha * (yct_W * yc)));
}

/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
 * excluding additive constants.
 *
 * @param phi Vector of parameters for spatial smoothing
 * @param n length of phi
 * @param node1
 * @param node2
 * @param w_ij weights corresponding to each pair of nodes
 * @param k number of groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param has_theta If the model contains an independent partial pooling term, phi for singletons can be zeroed out; otherwise, they require a standard normal prior. Both BYM and BYM2 have theta.
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, int n, int n_edges,
		      int[] node1, int[] node2, vector weight,
		      int k, int[] group_size, int[] group_idx,
		      int has_theta) {
  real lp;
  int pos=1;
  vector[n_edges] d = phi[node1] - phi[node2];
  lp = -0.5 * dot_product(weight, d .* d);
  if (has_theta) {
    for (j in 1:k) {
      /* sum to zero constraint, for each connected group; singletons zero out */
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
	/* its a singleton: independent */
	lp += std_normal_lpdf(phi[segment(group_idx, pos, group_size[j])]);
      }      
      pos += group_size[j];
    }
  }
  return lp;
}

/**
 * Combine local and global partial-pooling components into the convolved BYM2 term.
 *
 * @param phi local (spatially autocorrelated) component
 * @param theta global component
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param logit_rho proportion of convolution that is spatially autocorrelated, logit transformed
 * @param scale_factor scaling factor to put the ICAR prior on unit scale
 *
 * @return Log probability density of CAR prior up to additive constant
 */
vector convolve_bym2(vector phi, vector theta,
		      int n, int k,
		      int[] group_size, int[] group_idx,
		      real rho, vector scale_factor
		      ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
    convolution[ segment(group_idx, pos, group_size[j]) ] =
      sqrt(rho * inv(scale_factor[j])) * phi[ segment(group_idx, pos, group_size[j]) ] +
      sqrt(1 - rho) * theta[ segment(group_idx, pos, group_size[j]) ];
    pos += group_size[j];
  }
  return convolution;
}


