.onAttach <- function(libname, pkgname) {
    Sys.setenv("_SP_STARTUP_MESSAGE_"="none")
    packageStartupMessage("This is geostan version ", utils::packageVersion("geostan"))
}
