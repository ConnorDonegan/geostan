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
              vector A_w,
              array[] int A_v,
              array[] int A_u,
              vector D_inv,
              real log_det_D_inv,
              vector lambda,
              int n) {
  vector[n] z = y - mu;
  real ztDz = (z .* D_inv)' * z;
  real ztAz = z' * csr_matrix_times_vector(n, n, A_w, A_v, A_u, z);
  real ldet_ImrhoC = sum(log1m(rho * lambda));  
  return 0.5 * (
        -n * log( 2 * pi() )
        -2 * n * log(tau)
        + log_det_D_inv
        + ldet_ImrhoC
        - (1 / tau^2) * (ztDz - rho * ztAz));
}
