
context("SA indices")
test_that("SA indices produce the same results as spdep when expected", {
                                       
    skip_on_cran()
    skip_if_not_installed("spdep")
    skip_if_not_installed("spatialreg")
    library(spdep)
    library(spatialreg)
    d=10
    t=0.0001
    data(georgia)
    W <- shape2mat(georgia, style = "W")
    C <- shape2mat(georgia, style = "B")
    nb <- poly2nb(georgia)
    lw.W <- nb2listw(nb, style = "W")    
    lw.C <- nb2listw(nb, style = "B")
    x <- sim_sar(w=W, rho = 0.7)
    n <- length(x)
    ## test moran coefficient, binary matrix C
    mc.spdep <- round(as.numeric(moran.test(x, lw.C)$estimate["Moran I statistic"]), d)
    mc.geostan <- mc(x, C, digits = d)
    expect_equal(mc.spdep, mc.geostan, tolerance = t)
    ## test moran coefficient, row-standardized matrix W
    mc.spdep <- round(as.numeric(moran.test(x, lw.W)$estimate["Moran I statistic"]), digits = d)
    mc.geostan <- mc(x, W, digits = d)
    expect_equal(mc.spdep, mc.geostan)
    ## test APLE (estimate of SAR SA parameter)
    x <- as.numeric(scale(x))
    aple.spdep <- round(spdep::aple(x, lw.W), d)
    aple.geostan <- geostan::aple(x, W, digits = d)
    expect_equal(aple.spdep, aple.geostan)    
})

