/**
 * Log probability density of the simultaneous autoregressive (SAR) model (spatial error model)
 *
 * @param y Process to model
 * @param mu Mean vector
 * @param sigma Scale parameter
 * @param rho Spatial dependence parameter
 * @param W Sparse representation of W (its non-zero values)
 * @param W_v Column indices for values in W
 * @param W_u Row starting indices for values in W
 * @param lambda Eigenvalues of W
 * @param n Length of y
 *
 * @return Log probability density of SAR model up to additive constant
*/
  real  sar_normal_lpdf(vector y,
		      vector mu,
		      real sigma,
		      real rho,
		      vector W_w,
		      array[] int W_v,
		      array[] int W_u,
		      vector lambda,
		      int n) {
    vector[n] z = y - mu;
    real tau = 1 / sigma^2;
    vector[n] ImrhoWz = z - csr_matrix_times_vector(n, n, rho * W_w, W_v , W_u , z);
    real zVz = tau * dot_self(ImrhoWz);
    real ldet_V = 2 * sum(log1m(rho * lambda)) - 2 * n * log(sigma);
    return  0.5 * ( -n * log(2 * pi()) + ldet_V - zVz );		 
 }    
