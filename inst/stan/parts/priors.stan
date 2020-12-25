 /**
   * Regularized Horseshoe Prior (Piironen and Vehtari)
   * @return A vector of coefficientes
   */
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
  * intrinsic autoregressive prior (Morris et al 2019)
  * @return log probability density of IAR prior model up to an additive constant
  */
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }


/** Conditaional Autoregressive Prior (M. Joseph)
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param C_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param C_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*C*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, int[,] C_sparse, vector D_sparse, vector lambda, int n, int C_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_C; // phi' * C
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_C = rep_row_vector(0, n);
      for (i in 1:C_n) {
        phit_C[C_sparse[i, 1]] = phit_C[C_sparse[i, 1]] + phi[C_sparse[i, 2]];
        phit_C[C_sparse[i, 2]] = phit_C[C_sparse[i, 2]] + phi[C_sparse[i, 1]];
      }
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_C * phi)));
  }
