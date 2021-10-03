.onAttach <- function(libname, pkgname) {
    packageStartupMessage("This is geostan version ", utils::packageVersion("geostan"))
}
