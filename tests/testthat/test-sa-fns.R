## produces warnings without spatialreg package; but spatialreg is not need for geostan.
CRAN=TRUE
if (!CRAN) {
library(spdep)
context("SA indices")
test_that("SA indices produce the same results as spdep when expected", {
    d=4
    t=0.0001
    data(sentencing)
    W <- shape2mat(ohio, style = "W")
    C <- shape2mat(ohio, style = "B")
    nb <- poly2nb(ohio)
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

}
