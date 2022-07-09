## comment "skip..." and un-comment "devtools::load_all(..." to run test
skip("Run interactively to avoid crashing R session.")
##devtools::load_all("~/dev/geostan")
library(rstan)
library(testthat)
context("RStan conformity")
###########################################
## this line overwrites previous results ##
writeLines(paste0("# Results from tests/testthat/test-rstan-conformity \n# Run date: ",
                  Sys.Date()),
           "test-rstan-conformity-results")
###########################################

test_that("geostan model results match raw rstan results: gaussian glm", {
    georgia$y <- log(georgia$rate.male)
    fit1 <- stan_glm(y ~ ICE + insurance + college,
                     data = georgia,
                     prior = list(
                         intercept = normal(0, 1),
                         beta = normal(
                             location = c(0, 0, 0),
                             scale = c(10, 10, 10)
                         ),
                         sigma = student_t(df = 10, location = 0, scale = 1)
                         ),
                     family = gaussian(),
                     chains = 1,
                     iter= 6e3
                     )
    
    Sys.sleep(1)
    
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
    pars <- c("intercept", "beta", "sigma")
    b = apply( as.matrix(fit1.b, pars = pars), 2, mean)
    b = as.numeric(b)
    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)
    cat("\n\n# Context: geostan results equal Rstan results for normal GLM",
        file = "test-rstan-conformity-results",
        append = TRUE)
    cat(paste0("\n# ", pars, " geostan: ", round(a, 3), " rstan: ", round(b, 3)),
        file = "test-rstan-conformity-results",
        append = TRUE)
})

rm(list = ls() )

test_that("geostan model results match raw rstan results: Poisson rates", {
    fit2 <- stan_glm(deaths.male ~ offset(log(pop.at.risk.male)) + college,
                     data = georgia,
                     re = ~ GEOID,
                     prior = list(
                         intercept = normal(0, 10),
                         beta = normal(0, 10),
                         alpha_tau = student_t(df=10, location = 0, scale = 1)
                         ),                                     
                     family = poisson(),
                     chains = 4,
                     iter = 2e3
                     )
        
    Sys.sleep(1)
    
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
                   chains = 4,
                   iter = 2e3
                   )
    
    a = fit2$summary$mean
    a = as.numeric(a)
    pars <- c("intercept", "beta", "alpha_tau")
    b = apply( as.matrix(fit2.b, pars = pars), 2, mean)
    b = as.numeric(b)
    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)
    cat("\n\n# Context: geostan results equal Rstan results for log-linear Poisson model",
        file = "test-rstan-conformity-results",
        append = TRUE)
    cat(paste0("\n# ", pars, " geostan: ", round(a, 3), " rstan: ", round(b, 3)),
        file = "test-rstan-conformity-results",
        append = TRUE)

    f2 <- spatial(fit2)$mean
    f2.b <- apply(as.matrix(fit2.b, pars = "alpha_re"), 2, mean)
    f2.b <- as.numeric(f2.b)
    for (i in 1:length(f2)) expect_equal(f2[i], f2.b[i], tol = 0.01)
})

rm(list = ls())

test_that("geostan model results match raw rstan results: Binomial rates", {
    
    fit3 <- stan_glm(cbind(deaths.male, pop.at.risk.male - deaths.male) ~ college + ICE,
                     re = ~ GEOID,
                     data = georgia,
                     family = binomial(),
                     prior = list(intercept = normal(0, 10),
                                  beta = normal(location = c(0, 0),
                                                    scale = c(10, 10)),
                                  alpha_tau = student_t(10, 0, 1)
                                  ),
                     centerx = TRUE
                     )
        
    Sys.sleep(4)
    
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
                   chains = 4,
                   iter = 2e3
                   )
    a = fit3$summary$mean
    a = as.numeric(a)
    pars <- c("intercept", "beta[1]", "beta[2]", "alpha_tau")
    b = apply( as.matrix(fit3.b, pars = pars), 2, mean)
    b = as.numeric(b)
    for (i in 1:5) expect_equal(a[1], b[1], tol = 0.01)
    cat("\n\n# Context: geostan results equal Rstan results for binomial model",
        file = "test-rstan-conformity-results",
        append = TRUE)
    cat(paste0("\n# ", pars, " geostan: ", round(a, 3), " rstan: ", round(b, 3)),
        file = "test-rstan-conformity-results",
        append = TRUE)
    
    f3 <- spatial(fit3)$mean
    f3.b <- apply(as.matrix(fit3.b, pars = "alpha_re"), 2, mean)
    f3.b <- as.numeric(f3.b)
    for (i in 1:length(f3)) expect_equal(f3[i], f3.b[i], tol = 0.01)
                  
})


rm(list=ls())
test_that("geostan model results match raw rstan results: non-spatial ME model", { 
    me_prior <- list(location = normal(location = 0,
                                       scale = 0.5),
                 scale = student_t(df = 10,
                                    location = 0,
                                    scale = 1)
                 )
    ME <- prep_me_data(se = data.frame(ICE = georgia$ICE.se),
                       prior = me_prior) 
    fit4.a <- stan_glm(log(rate.male) ~ ICE,
                     data = georgia,
                     ME = ME,
                     prior = list(intercept = normal(0, 10),
                                  beta = normal(location = 0, scale = 2),
                                  sigma = student_t(10, 0, 1)
                                  ),
                     chains = 4,
                     iter = 2e3
                     )
    
    Sys.sleep(1)
        
mod4 <- "
data {
    int n;
    vector[n] y;
    vector[n] z;
    vector[n] se;
}
parameters {
    vector[n] x;
    real<lower=0> df_x;
    real mu_x;
    real<lower=0> sigma_x;
    real intercept;
    real beta;
    real<lower=0> sigma;
}
transformed parameters {
}

model {
    target += normal_lpdf(z | x, se);
    target += student_t_lpdf(x | df_x, mu_x, sigma_x);
    target += normal_lpdf(y | intercept + x * beta, sigma);
    target += normal_lpdf(intercept | 0, 10);
    target += normal_lpdf(beta | 0, 2);
    target += student_t_lpdf(sigma | 10, 0, 1);
    mu_x ~ normal(0, 0.5);
    df_x ~ gamma(3, 0.2);
    sigma_x ~ student_t(10, 0, 1);
}

"
    m4 <- rstan::stan_model(model_code = mod4)
    fit4.b <- rstan::sampling(m4,
                   data = list(n = nrow(georgia),
                               y = log(georgia$rate.male),
                               z = georgia$ICE,
                               se = georgia$ICE.se
                               ),
                   chains = 4,
                   iter = 2e3
                   )
                                        # nu_x_true
    pars <- c("intercept", "beta", "sigma", "nu_x_true", "mu_x_true", "sigma_x_true")
    a <- apply( as.matrix(fit4.a, pars = pars), 2, mean)
    a = as.numeric(a)
    b = apply( as.matrix(fit4.b, pars = c("intercept", "beta", "sigma", "df_x", "mu_x", "sigma_x")), 2, mean)
    b = as.numeric(b)
    for (i in 1:length(a)) expect_equal(a[1], b[1], tol = 0.01)
    cat("\n\n# Context: geostan results equal Rstan results for ME model",
        file = "test-rstan-conformity-results",
        append = TRUE)
    cat(paste0("\n# ", pars, " geostan: ", round(a, 3), " rstan: ", round(b, 3)),
        file = "test-rstan-conformity-results",
        append = TRUE)
})
