# R/zzz.R

# List the specific callable routines (your BIND(C, NAME) strings)
# This is the most minimal and reliable registration list.
.FortranRoutines <- c("twcdf", "twpdf")

# Explicitly register the routines when the package loads
.onLoad <- function(libname, pkgname) {
  # This uses R's internal mechanism to find and register the compiled Fortran routines.
  # We use the full, safer set of names here to guarantee linkage.
  
  # Note: The DLL name is just the package name (e.g., "tweedie")
  # We use the official R function for native routine registration.
  
  # Register the main function (or functions) called from R:
  # This forces R to check the DLL for these specific entry points.
  if (requireNamespace("utils", quietly = TRUE)) {
    utils::assignNativeRoutines(
      getLoadedDLLs()[[pkgname]], # Finds the loaded DLL object for your package
      routines = .FortranRoutines
    )
  }
}

.onUnload <- function(libpath) {
  # Unload the dynamic library when the package is detached
  library.dynam.unload("tweedie", libpath)
}