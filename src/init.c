#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Forward declaration for the ONLY top-level routine called by R. */
// This must match the BIND(C, NAME='twcdf') signature (9 arguments).
extern void twcdf(
    int *N, 
    double *power, 
    double *phi, 
    double *y, 
    double *mu, 
    double *cdf_out, 
    int *exitstatus, 
    double *relerr, 
    int *its
); 

/* Registering ONLY the top-level R entry point. */
// Internal routines (DFbigp, DFsmallp, rtnewton, etc.) will be resolved dynamically.
static const R_FortranMethodDef FortranEntries[] = {
  // Top-level R entry point (9 arguments)
  {"twcdf", (DL_FUNC) &twcdf, 9}, 
  {NULL, NULL, 0}
};

void R_init_tweedie(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  // CRITICAL: Ensure dynamic symbol resolution is possible for internal calls.
  R_useDynamicSymbols(dll, FALSE); 
}
