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
 * Log probability density of the simultaneous autoregressive (SAR) model (spatial error model)
 *
 * @param y Process to model
 * @param mu Mean vector
 * @param sigma Scale parameter
 * @param rho Spatial dependence parameter
 * @param ImW Sparse representation of (I - W): non-zero values only
 * @param ImW_v Column indices for values in ImW
 * @param ImW_u Row starting indices for values in ImW
 * @param Widx Indices for the off-diagonal elements in ImC
 * @param lambda Eigenvalues of W
 * @param n Length of y
 *
 * @return Log probability density of SAR model up to additive constant
*/
  real  sar_normal_lpdf(vector y,
		      vector mu,
		      real sigma,
		      real rho,
		      vector ImW,
		      int[] ImW_v,
		      int[] ImW_u,
		      int[] Widx,
		      vector lambda,
		      int n) {
  vector[n] z = y - mu;
  real tau = 1 / sigma^2;
  vector[num_elements(ImW)] ImrhoW = ImW;  // (I - rho W)
  vector[n] ImrhoWz;                       // (I - rho * W) * z
  real zVz;
  real ldet_V = 2 * sum(log1m(rho * lambda)) - 2 * n * log(sigma);
  ImrhoW[Widx] = rho * ImW[Widx];
  ImrhoWz = csr_matrix_times_vector(n, n, ImrhoW, ImW_v , ImW_u , z);
  zVz = tau * dot_self(ImrhoWz);
  return  0.5 * (-n * log(2 * pi()) + ldet_V - zVz);		 
 }    


/**
 * Log probability density of the conditional autoregressive (CAR) model
 *
 * @param y Process to model
 * @param mu Mean vector
 * @param tau Scale parameter
 * @param rho Spatial dependence parameter
 * @param ImC Sparse representation of (I - C): non-zero values only
 * @param ImC_v Column indices for values in ImC
 * @param ImC_u Row starting indices for values in ImC
 * @param Cidx Indices for the off-diagonal elements in ImC
 * @param D_inv Diagonal elements from the inverse of Delta, where M = Delta * tau^2 is a diagonal matrix containing the conditional variances.
 * @param log_det_D_inv Log determinant of Delta inverse
 * @param lambda Eigenvalues of C (or of the symmetric, scaled matrix Delta^{-1/2}*C*Delta^{1/2}).
 * @param n Length of y
 *
 * @return Log probability density of CAR model up to additive constant
*/
real car_normal_lpdf(vector y, vector mu,
		     real tau, real rho,
		     vector ImC, int[] ImC_v, int[] ImC_u, int[] Cidx,
		     vector D_inv, real log_det_D_inv, vector lambda,
		     int n) {
  vector[n] z = y - mu;
  vector[num_elements(ImC)] ImrhoC = ImC; // (I - rho C)
  vector[n] zMinv = (1 / tau^2) * z .* D_inv; // z' * M-1
  vector[n] ImrhoCz; // (I - rho * C) * z
  vector[n] ldet_ImrhoC;
  ImrhoC[Cidx] = rho * ImC[Cidx];
  ImrhoCz = csr_matrix_times_vector(n, n, ImrhoC, ImC_v, ImC_u, z);
  for (i in 1:n) ldet_ImrhoC[i] = log1m(rho * lambda[i]);
  return 0.5 * (
		-n * log( 2 * pi() )
		- 2 * n * log(tau)
		+ log_det_D_inv
		+ sum(ldet_ImrhoC)
		- dot_product(zMinv, ImrhoCz)
		);
}

/**
 * Log probability density of the conditional autoregressive (CAR) model: WCAR specifications only
 *
 * @param y Process to model
 * @param mu Mean vector
 * @param tau Scale parameter
 * @param rho Spatial dependence parameter
 * @param A_w Sparse representation of the symmetric connectivity matrix, A
 * @param A_v Column indices for values in A_w
 * @param A_u Row starting indices for values in A_u
 * @param D_inv The row sums of A; i.e., the diagonal elements from the inverse of Delta, where M = Delta * tau^2 is a diagonal matrix containing the conditional variances.
 * @param log_det_D_inv Log determinant of Delta inverse.
 * @param lambda Eigenvalues of C (or of the symmetric, scaled matrix Delta^{-1/2}*C*Delta^{1/2}); for the WCAR specification, C is the row-standardized version of W.
 * @param n Length of y
 *
 * @return Log probability density of CAR prior up to additive constant
 */
real wcar_normal_lpdf(vector y, vector mu,
		      real tau, real rho,
		      vector A_w, int[] A_v, int[] A_u,
		      vector D_inv, real log_det_D_inv,
		      vector lambda,
		      int n) {
  vector[n] z = y - mu;
  real ztDz; // z transpose * D * z
  real ztAz; // z transpose * A * z
  vector[n] ldet_ImrhoC;
  ztDz = (z .* D_inv)' * z;
  ztAz = z' * csr_matrix_times_vector(n, n, A_w, A_v, A_u, z);
  for (i in 1:n) ldet_ImrhoC[i] = log1m(rho * lambda[i]);
  return 0.5 * (
		-n * log( 2 * pi() )
		-2 * n * log(tau)
		+ log_det_D_inv
		+ sum(ldet_ImrhoC)
		- (1 / tau^2) * (ztDz - rho * ztAz));
}

/** 
 * Conditional Autoregressive Model
 */
real auto_normal_lpdf(vector y, vector mu,
		      real tau, real rho,
		      vector Ax_w, int[] Ax_v, int[] Ax_u,
		      int[] Cidx,
		      vector D_inv, real log_det_D_inv,
		      vector lambda,
		      int n, int WCAR) {
  if (WCAR) {
    return wcar_normal_lpdf(y | mu, tau, rho, Ax_w, Ax_v, Ax_u, D_inv, log_det_D_inv, lambda, n);
  } else {
    return car_normal_lpdf(y | mu, tau, rho, Ax_w, Ax_v, Ax_u, Cidx, D_inv, log_det_D_inv, lambda, n);
      }
}

/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
 * excluding additive constants. 
 *
 * @param phi Vector of parameters for spatial smoothing (on unit scale)
 * @param spatial_scale Scale parameter for the ICAR model
 * @param node1 
 * @param node2
 * @param k number of groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 * @param has_theta If the model contains an independent partial pooling term, phi for singletons can be zeroed out; otherwise, they require a standard normal prior. Both BYM and BYM2 have theta.
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, real spatial_scale,
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
    /* does not have theta */
    for (j in 1:k) {
      if (group_size[j] > 1) {
    /* same as above for non-singletons: sum to zero constraint */
    lp += normal_lpdf(sum(phi[segment(group_idx, pos, group_size[j])]) | 0, 0.001 * group_size[j]);
      } else {
    /* its a singleton: independent Gaussian prior on phi */
    lp += normal_lpdf(phi[ segment(group_idx, pos, group_size[j]) ] | 0, spatial_scale);
      }      
      pos += group_size[j];
    }
  }
  return lp;
}

/**
 * Create phi from phi_tilde, inv_sqrt_scale_factor, and spatial_scale. 
 *
 * @param phi_tilde local component (spatially autocorrelated) 
 * @param phi_scale scale parameter for phi
 * @param rho proportion spatial (for ICAR and BYM models, this alwasy equals 1; for BYM2, it is a model parameter)
 * @param inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 *                              transformed from 1/scale^2 --> scale. Or, a vector of ones.
 * @param n number of spatial units
 * @param k number of connected groups
 * @param group_size number of observational units in each group
 * @param group_idx index of observations in order of their group membership
 *
 * @return phi vector of spatially autocorrelated coefficients
 */
vector make_phi(vector phi_tilde, real phi_scale,
		real rho,
		vector inv_sqrt_scale_factor,
		int n, int k,
		int[] group_size, int[] group_idx
              ) {
  vector[n] phi;
  int pos=1;
  for (j in 1:k) {
    phi[ segment(group_idx, pos, group_size[j]) ] = phi_scale * sqrt(rho) * inv_sqrt_scale_factor[j] * phi_tilde[ segment(group_idx, pos, group_size[j]) ];
    pos += group_size[j];
  }
  return phi;
}


/**
 * Combine local and global partial-pooling components into the convolved BYM term.
 *
 * @param phi spatially autocorrelated component (not phi_tilde!)
 * @param theta global component (not theta_tilde!)
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
 * @param inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 *                              transformed from 1/scale^2 --> scale. Or, a vector of ones.
 *
 * @return BYM2 convolution vector
 */
vector convolve_bym2(vector phi_tilde, vector theta_tilde,
          real spatial_scale,
		      int n, int k,
		      int[] group_size, int[] group_idx,
		      real rho, vector inv_sqrt_scale_factor
		      ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:k) {
    if (group_size[j] == 1) {
        convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * theta_tilde[ segment(group_idx, pos, group_size[j]) ];
    } else {
    convolution[ segment(group_idx, pos, group_size[j]) ] = spatial_scale * (
     sqrt(rho) * inv_sqrt_scale_factor[j] * phi_tilde[ segment(group_idx, pos, group_size[j]) ] +
     sqrt(1 - rho) * theta_tilde[ segment(group_idx, pos, group_size[j]) ]
      );
  }
  pos += group_size[j];
  }
  return convolution;
}


