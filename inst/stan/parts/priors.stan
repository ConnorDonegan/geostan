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

