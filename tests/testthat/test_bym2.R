
iter=10
refresh=0
source("helpers.R")

##devtools::load_all("~/dev/geostan")

context("stan_bym2")
test_that("BYM2 with offset model works", {
#    library(INLA)
    data(sentencing)
    n <- nrow(sentencing)
    C <- shape2mat(sentencing)
    ME <- list(offset = rep(10, n))
    ## nbs <- edges(C)
    ## N <- nrow(C)
    ## adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
    ## Q=  Diagonal(N, rowSums(adj.matrix)) - adj.matrix
    ## Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
    ## Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))
    ## scaling_factor = exp(mean(log(diag(Q_inv))))
    scaling_factor = 0.71
    SW(
        fit <- stan_bym2(sents ~ offset(expected_sents),
                    data = sentencing,
                  #  ME = ME,
                    C = C,
                    scaleFactor = scaling_factor, 
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("BYM2 works with covariate ME", {
    data(sentencing)
    n <- nrow(sentencing)
    sentencing$x <- rnorm(n = n, sd = 1)
    C <- shape2mat(sentencing)
    ME <- list(offset = rep(10, n),
               ME = data.frame(x = rep(.1, n)))
    scaling_factor = 0.71
    SW(
        fit <- stan_bym2(sents ~ offset(expected_sents) + x,
                    data = sentencing,
                    slx = ~ x,
                    ME = ME,
                    C = C,
                    scaleFactor = scaling_factor, 
                    chains = 1,
                    family = poisson(),
                    iter = iter,
                    refresh = refresh)
    )
    expect_geostan(fit)
})

test_that("BYM2 binomial works and accepts covariate ME", {
    data(ohio)
    C <- shape2mat(ohio)        
    n <- nrow(ohio)
    ME <- list(ME = data.frame(unemployment = rep(0.25, n),
                               historic_gop = rep(1, n),
                               college_educated = rep(1, n)),
               bounded = c(1, 0, 1))
    scaling_factor = 0.40
    SW(
        fit <- stan_bym2(cbind(trump_2016, total_2016 - trump_2016) ~ log(population) + college_educated + unemployment + historic_gop,
                         data = ohio,
                         family = binomial(),
                         C = C,
                         ME = ME,
                         scaleFactor = scaling_factor,
                         iter = iter,
                         refresh = refresh,
                         chains = 1)
                         )
                         
    expect_geostan(fit)
})

