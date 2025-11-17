#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h> 

// (Routines are extern void-declared here, as you already have them)

/* --- C and Fortran Routine Registration Tables here --- */
// (CEntries and FortranEntries tables go here, as you already have them)


// --- The R-MANDATORY initialization entry point ---
// The function must be exported and registered. 
// We rely on the R build system to handle the visibility, not non-portable attributes.
void R_init_tweedie(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL); 
  R_useDynamicSymbols(dll, FALSE); 
}