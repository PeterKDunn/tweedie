#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> // <-- New required header

/* --- Forward declarations for BIND(C) Fortran Subroutines --- */
extern void twcomputation(int *N, double *p, double *phi, double *y, double *mu, int *verbose, int *pdf, double *funvalue, int *exitstatus, double *relerr, int *its);

/* --- C Routine Registration Table (R_CMethodDef) --- */
static const R_CMethodDef CEntries[] = {
  {"twcomputation", (DL_FUNC) twcomputation, 11, {R_TYPE_INT, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_INT, R_TYPE_INT, R_TYPE_DBL, R_TYPE_INT, R_TYPE_DBL, R_TYPE_INT}},
  {NULL, NULL, 0, {0}}
};

/* --- Fortran Routine Registration Table (R_FortranMethodDef) --- */
static const R_FortranMethodDef FortranEntries[] = {
  {NULL, NULL, 0}
};


/* --- The main initialization entry point --- */
// Force the function to be exported (visible) on Linux
void R_init_tweedie(DllInfo *dll) attribute_visible; 

void R_init_tweedie(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL); 
  R_useDynamicSymbols(dll, FALSE); 
}