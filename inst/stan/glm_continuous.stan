data { 
#include parts/glm_data.stan
  vector[n] y;
}

transformed data { 
#include parts/QR.stan
}

parameters {
#include parts/glm_parameters.stan
}

transformed parameters {
#include parts/glm_trans_params.stan
}

model {
#include parts/glm_model.stan
#include parts/y_continuous_model.stan
}

generated quantities {
  vector[n] log_lik;
  vector[n] yrep;
  vector[n] residual;
  vector[n] fitted;
  vector[n_ids] alpha_re;
  if (has_re) {
    for (i in 1:n_ids) {
      alpha_re[i] = alpha_tau[has_re] * alpha_re_tilde[i];
    }
  }
  for (i in 1:n) {
   fitted[i] = f[i];
   residual[i] = y[i] - fitted[i];
   if (is_student) {
      log_lik[i] = student_t_lpdf(y[i] | nu[1], fitted[i], sigma[1]);
      yrep[i] = student_t_rng(nu[1], fitted[i], sigma[1]);     
   } else {
      log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma[1]);
      yrep[i] = normal_rng(fitted[i], sigma[1]);     
   }
 }
}

