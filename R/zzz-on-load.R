.onAttach <- function(libname, pkgname) {
    packageStartupMessage("This is a development version of geostan version ", utils::packageVersion("geostan"))
}
