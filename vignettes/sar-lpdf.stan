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
              array[] int ImW_v,
              array[] int ImW_u,
              array[] int Widx,
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
 
 
