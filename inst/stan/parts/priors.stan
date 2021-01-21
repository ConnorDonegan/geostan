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

/**
 * Return the log probability of a conditional autoregressive (CAR) model,
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
