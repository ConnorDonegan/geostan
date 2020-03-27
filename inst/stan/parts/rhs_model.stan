  z ~ std_normal();
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5, 0.5);// .5 * nu_local, .5 * nu_local, nu_local = 1
  aux1_global ~ std_normal();
  aux2_global ~ inv_gamma(0.5, 0.5); // .5 * nu_local, .5 * nu_global, both = 1
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);

