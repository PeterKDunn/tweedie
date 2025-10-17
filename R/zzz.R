# R/zzz.R

.onLoad <- function(libname, pkgname) {
}

.onUnload <- function(libpath) {
  # library.dynam() manually in .onLoad.
  library.dynam.unload("tweedie", libpath)
}