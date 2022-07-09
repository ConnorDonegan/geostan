
context("SA indices")
test_that("SA indices produce the same results as spdep when expected", {
    skip_on_cran()
#    skip_if_not_installed("spatialreg")
    d=10
    data(georgia)
    W <- shape2mat(georgia, style = "W")
    C <- shape2mat(georgia, style = "B")
    nb <- spdep::poly2nb(georgia)
    lw.W <- spdep::nb2listw(nb, style = "W")    
    lw.C <- spdep::nb2listw(nb, style = "B")
    x <- sim_sar(w=W, rho = 0.7)
    n <- length(x)
    ## Moran coefficient, binary matrix C
    mc_spdep <- spdep::moran.test(x, lw.C)$estimate["Moran I statistic"]
    mc_spdep <- as.numeric(mc_spdep)
    mc_geostan <- mc(x, C, digits = d)
    expect_equal(mc_geostan, mc_spdep)
    ## Moran coefficient, row-standardized matrix W
    mc_spdep <- spdep::moran.test(x, lw.W)$estimate["Moran I statistic"]
    mc_spdep <- as.numeric(mc_spdep)
    mc_geostan <- mc(x, W, digits = d)
    expect_equal(mc_geostan, mc_spdep)
    ## Local Moran's I: is not calculated in the same manner
    geo_lisa <- geostan::lisa(x, W, type = FALSE, digits = d)
    spdep_lisa <- as.numeric(spdep::localmoran(x, lw.W)[,"Ii"])
    m2 <- sum(scale(x, center=T, scale = T)^2)/length(x)
    expect_equal(geo_lisa/m2, spdep_lisa)
    ## APLE (estimate of SAR SA parameter)
    x <- as.numeric(scale(x))
    aple_spdep <- spdep::aple(x, lw.W)
    aple_geostan <- geostan::aple(x, W, digits = d)
    expect_equal(aple_spdep, aple_geostan)    
    ## Geary Ratio
    GR_geostan <- gr(x, W, digits = d)
    GR_spdep <- spdep::geary.test(x, lw.W)
    expect_equal(GR_geostan, as.numeric(GR_spdep$estimate["Geary C statistic"]))    
    ## Local Geary
    LG_geostan <- lg(x, W, digits = d)
    LG_spdep <- spdep::localC(x, lw.W)
    expect_equal(LG_geostan, LG_spdep)
})

