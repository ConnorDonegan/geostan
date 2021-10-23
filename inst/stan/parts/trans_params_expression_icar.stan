  if (type == 1) {
    phi = make_phi(phi_tilde, spatial_scale[1], 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += phi;
  }
  if (type == 2) {
    theta = theta_tilde * theta_scale[1];
    phi = make_phi(phi_tilde, spatial_scale[1], 1, inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym(phi, theta, n, k, group_size, group_idx);
  }  
  if (type == 3) {
    theta = spatial_scale[1] * sqrt(1 - rho[1]) * theta_tilde;
    phi = make_phi(phi_tilde, spatial_scale[1], rho[1], inv_sqrt_scale_factor, n, k, group_size, group_idx);
    if (m) phi += A * alpha_phi;    
    fitted += convolve_bym2(phi_tilde, theta_tilde, spatial_scale[1], n, k, group_size, group_idx, rho[1], inv_sqrt_scale_factor);    
  }
