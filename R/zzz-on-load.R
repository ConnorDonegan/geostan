.onAttach <- function(libname, pkgname) {
    packageStartupMessage("This is geostan version ", utils::packageVersion("geostan"))
    packageStartupMessage("This package is still under development and has not been released for general use.")
}
