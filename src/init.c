#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Forward declarations for the Fortran subroutines (11 arguments) */
void twcdf(void*, void*, void*, void*, void*, void*, void*, 
              void*, void*, void*);
void twpdf(void*, void*, void*, void*, void*, void*, void*, 
              void*, void*, void*, void*);

/* List of Fortran routines callable from R */
static const R_FortranMethodDef FortranEntries[] = {
  {"twcdf", (DL_FUNC) twcdf, 9}, 
  {"twpdf", (DL_FUNC) twpdf, 10},
  {NULL, NULL, 0}
};

/* The main initialization entry point R looks for */
void R_init_tweedie(DllInfo *dll)
{
  // The 5 expected lists are: C, Call, Fortran, External.
  // DllInfo is the first argument, so 4 more are expected.
  // The final call MUST have only 5 lists (1 DllInfo + 4 routine lists).
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL); 
}