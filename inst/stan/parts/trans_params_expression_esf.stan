  if (dev) {
    if (has_sigma) {
      error_scale[1] = sigma[1];
    } else {
      error_scale[1] = 1;
    }
    beta_ev = rhs_prior(dev, z, aux1_global[1], aux2_global[1], aux1_local, aux2_local, caux[1], global_scale, slab_scale, error_scale[1]);
    esf = EV * beta_ev;
    fitted += esf;
  }
