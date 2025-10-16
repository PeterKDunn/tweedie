# R/zzz.R

# (Keep your .FortranRoutines list)
.FortranRoutines <- c(
  "twpdf",    
  "twcdf"
)

.onLoad <- function(libname, pkgname) {
  # 1. Explicitly load the dynamic library (guarantees the file is mapped to memory)
  # This call returns the DLLInfo object.
  # We assign it to a local variable to confirm loading.
  DLL <- library.dynam(pkgname, pkgname, libname)
  
  # 2. Check that the DLL object was successfully loaded
  if (is.null(DLL)) {
    warning("Failed to load DLL for package ", pkgname)
  }
  
  # Note: The presence of R_init_tweedie.f90 and NAMESPACE .registration=TRUE 
  # handles the routine registration once the DLL is loaded successfully above.
}

.onUnload <- function(libpath) {
  # Standard cleanup
  library.dynam.unload("tweedie", libpath)
}