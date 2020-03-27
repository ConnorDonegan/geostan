  if (is_student) {
    nu[1] ~ gamma(t_nu_prior[1], t_nu_prior[2]);
    y ~ student_t(nu[1], f, sigma[1]);
  } else {
    y ~ normal(f, sigma[1]);
  }

