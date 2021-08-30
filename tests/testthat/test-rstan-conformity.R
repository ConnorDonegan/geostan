####devtools::load_all("~/dev/geostan")
library(rstan)


context("RStan conformity")
skip("Run these interactively to avoid crashing R session.")
test_that("geostan model results match raw rstan results: gaussian glm", {

    georgia$y <- log(georgia$rate.male)
    fit1 <- stan_glm(y ~ ICE + insurance + college,
                     data = georgia,
                     prior = list(
                         intercept = c(0, 1),
                         beta = data.frame(
                             location = c(0, 0, 0),
                             scale = c(10, 10, 10)
                         ),
                         sigma = c(df = 10, location = 0, scale = 1)
                     ),                 
                     family = gaussian(),
                     chains = 1,
                     iter= 6e3
                     )
    
    mod1 <- "
data {
    int n;
    int k;
    vector[n] y;
    matrix[n, k] x;
}
parameters {
    real intercept;
    vector[k] beta;
    real<lower=0> sigma;
}
transformed parameters {
  vector[n] fitted = intercept + x * beta;
}
model {
    target += normal_lpdf(y | fitted, sigma);
    target += student_t_lpdf(sigma | 10, 0, 1);
    target += normal_lpdf(intercept | 0, 1);
    target += normal_lpdf(beta | 0, 10);
}
generated quantities {
 vector[n] log_lik;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i] | fitted[i], sigma);
}

"
    m1 <- rstan::stan_model(model_code = mod1)


    fit1.b <- stan(model_code = mod1,
                   data = list(n = nrow(georgia),
                               k = 3,
                               y = georgia$y,
                               x = model.matrix(~ 0 + ICE + insurance + college, georgia)
                               ),
                   chains = 1,
                   iter = 6e3
                   )


    a = fit1$summary$mean
    a = as.numeric(a)

    b = apply( as.matrix(fit1.b, pars = c("intercept", "beta", "sigma")), 2, mean)
    b = as.numeric(b)

    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)
})

rm(list = ls() )

test_that("geostan model results match raw rstan results: Poisson rates", {
    
    fit2 <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + college,
                     data = georgia,
                     re = ~ GEOID,
                     priot = list(
                         intercept = c(0, 10),
                         beta = c(0, 10),
                         alpha_tau = c(df=10, location = 0, scale = 1)
                         ),                                     
                     family = poisson(),
                     chains = 1,
                     iter= 6e3
                     )
    
mod2 <- "
data {
    int n;
    int k;
    int y[n];
    matrix[n, k] x;
    vector[n] log_at_risk;
}
parameters {
    real intercept;
    vector[k] beta;
    vector[n] alpha_re;
    real<lower=0> alpha_tau;
}
transformed parameters {
  vector[n] fitted = exp(log_at_risk + intercept + x * beta + alpha_re);
}
model {
    target += poisson_lpmf(y | fitted);
    target += normal_lpdf(intercept | 0, 10);
    target += normal_lpdf(beta | 0, 10);
    target += normal_lpdf(alpha_re | 0, alpha_tau);
    target += student_t_lpdf(alpha_tau | 10, 0, 1);
}
generated quantities {
 vector[n] log_lik;
  for (i in 1:n) log_lik[i] = poisson_lpmf(y[i] | fitted[i]);
}

"

    m2 <- rstan::stan_model(model_code = mod2)

    fit2.b <- rstan::sampling(m2,
                   data = list(n = nrow(georgia),
                               k = 1,
                               y = georgia$deaths.male,
                               x = model.matrix(~ 0 + college, georgia),
                               log_at_risk = log(georgia$pop.at.risk.male)
                               ),
                   chains = 1,
                   iter = 6e3
                   )


    a = fit2$summary$mean
    a = as.numeric(a)

    b = apply( as.matrix(fit2.b, pars = c("intercept", "beta", "alpha_tau")), 2, mean)
    b = as.numeric(b)

    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)

    f2 <- spatial(fit2)$mean
    f2.b <- apply(as.matrix(fit2.b, pars = "alpha_re"), 2, mean)
    f2.b <- as.numeric(f2.b)

    for (i in 1:length(f2)) expect_equal(f2[i], f2.b[i], tol = 0.01)
                  
})

rm(list = ls() )
test_that("geostan model results match raw rstan results: Binomial rates", {


    #### slower than the stan model. I bet that the RE term should be centered on alpha, not zero.
    fit3 <- stan_glm(cbind(deaths.male, pop.at.risk.male - deaths.male) ~ college + ICE,
                     re = ~ GEOID,
                     data = georgia,
                     family = binomial(),
                     prior = list(intercept = c(0, 10),
                                  beta = data.frame(location = c(0, 0),
                                                    scale = c(10, 10)),
                                  alpha_tau = c(10, 0, 1)
                                  ),
                     centerx = TRUE,
                     chains = 1,
                     iter= 8e3
                     )
    
mod3 <- "
data {
    int n;
    int k;
    int y[n];
    int trials[n];
    matrix[n, k] x;
}
parameters {
    real intercept;
    vector[k] beta;
    vector[n] alpha_re;
    real<lower=0> alpha_tau;
}
transformed parameters {
  vector[n] fitted = inv_logit(intercept + x * beta + alpha_re);
}
model {
    target += binomial_lpmf(y | trials, fitted);
    target += normal_lpdf(intercept | 0, 10);
    target += normal_lpdf(beta | 0, 10);
    target += normal_lpdf(alpha_re | 0, alpha_tau);
    target += student_t_lpdf(alpha_tau | 10, 0, 1);
}
generated quantities {
 vector[n] log_lik;
  for (i in 1:n) log_lik[i] = binomial_lpmf(y[i] | trials[i], fitted[i]);
}

"

    m3 <- rstan::stan_model(model_code = mod3)

    fit3.b <- rstan::sampling(m3,
                   data = list(n = nrow(georgia),
                               k = 2,
                               y = georgia$deaths.male,
                               trials = georgia$pop.at.risk.male,
                               x = scale(model.matrix(~ 0 + college + ICE, georgia), center = T, scale = F)
                               ),
                   chains = 1,
                   iter = 8e3
                   )


    a = fit3$summary$mean
    a = as.numeric(a)
    b = apply( as.matrix(fit3.b, pars = c("intercept", "beta", "alpha_tau")), 2, mean)
    b = as.numeric(b)
    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)

    f3 <- spatial(fit3)$mean
    f3.b <- apply(as.matrix(fit3.b, pars = "alpha_re"), 2, mean)
    f3.b <- as.numeric(f3.b)
    for (i in 1:length(f3)) expect_equal(f3[i], f3.b[i], tol = 0.01)
                  
})
