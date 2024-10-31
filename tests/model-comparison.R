##
## spatial econometric models and model comparison with DIC, WAIC, and model probabilities. N = 520
##

#devtools::load_all("~/dev/geostan")

library(geostan)

row = 20
col = 26
N <- row * col
sdl <- prep_sar_data2(row = row, col = col, quiet = TRUE)
w <- sdl$W
x <- sim_sar(w = w, rho = .5)
mu = (.5 * x - .25 * w %*% w)[,1]
y <- sim_sar(w = w, rho = .5, mu = mu)
df <- data.frame(y=y, x=x)

iter = 175

sdlm <- stan_sar(y ~ x,
                data = df,
                type = "SDLM",
                sar_parts = sdl,
                chains = 1,
                iter = iter,
                quiet = TRUE,
                keep_all = TRUE,
                slim = F)

sdem <- stan_sar(y ~ x,
                 data = df,
                 type = "SDEM",
                 sar_parts = sdl,
                 chains = 1,
                 iter = iter,
                 quiet = TRUE,
                keep_all = TRUE,                 
                slim = F) 

print(rbind(dic(sdlm), dic(sdem)))
print(rbind(waic(sdlm), waic(sdem)))


  #  library(bridgesampling)
  #  sdlm <- bridge_sampler(sdlm$stanfit)
  #  sdem <- bridge_sampler(sdem$stanfit)
  #  cat("post_prob (SDLM, SDEM):\n", post_prob(sdlm, sdem), "\n")

