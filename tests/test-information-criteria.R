devtools::load_all("~/dev/geostan")
row = 20
col = 26
N <- row * col
sdl <- prep_sar_data2(row = row, col = col, quiet = TRUE)
x <- rnorm(n = N)
y <- .75 *x + rnorm(n = N, sd = .5)
df <- data.frame(y=y, x=x)

sdlm <- stan_sar(y ~ x,
                data = df,
                type = "SDLM",
                sar_parts = sdl,
                chains = 1,
                iter = 500,
                quiet = TRUE,
                slim = F)


sdem <- stan_sar(y ~ x,
                 data = df,
                 type = "SDEM",
                 sar_parts = sdl,
                 chains = 1,
                 iter = 500,
                 quiet = TRUE,
                 slim = F)

print(rbind(dic(sdlm), dic(sdem)))
print(rbind(waic(sdlm), waic(sdem)))
