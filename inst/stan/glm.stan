data { 
#include parts/data.stan
}

transformed data {
#include parts/trans_data.stan
}

parameters {
#include parts/params.stan
}

transformed parameters {
#include parts/trans_params_declaration.stan
#include parts/trans_params_expression.stan
}

model {
#include parts/model.stan
}

generated quantities {
#include parts/gen_quants_declaration.stan
  for (i in 1:n) {
#include parts/gen_quants_expression_in_loop.stan
 }
}

