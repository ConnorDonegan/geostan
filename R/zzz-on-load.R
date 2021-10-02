.onAttach <- function(libname, pkgname) {
    packageStartupMessage("This is geostan version ", utils::packageVersion("geostan"))
    packageStartupMessage("To report issues, or find help, visit https://github.com/ConnorDonegan/geostan/issues.")
}
