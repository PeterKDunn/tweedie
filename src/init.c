#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* --- Forward declarations for BIND(C) Fortran Subroutines --- */
// twcomputation has 11 arguments
extern void twcomputation(int *N, double *p, double *phi, double *y, double *mu, int *verbose, int *pdf, double *funvalue, int *exitstatus, double *relerr, int *its);
/* --- C Routine Registration Table (R_CMethodDef) --- */
/* R will look here for BIND(C) symbols like "twcomputation" */
static const R_CMethodDef CEntries[] = {
  // R name    C pointer (The BIND(C) name) Arg Count   Argument Types
  {"twcomputation", (DL_FUNC) twcomputation, 10, {R_TYPE_INT, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_DBL, R_TYPE_INT, R_TYPE_DBL, R_TYPE_INT, R_TYPE_DBL, R_TYPE_INT}},
};


/* --- Fortran Routine Registration Table (R_FortranMethodDef) --- */
/* This table is empty now, as BIND(C) routines go in CEntries */
static const R_FortranMethodDef FortranEntries[] = {
  {NULL, NULL, 0, {0}} // EMPTY TABLE
};


/* --- The main initialization entry point R looks for (Single Definition) --- */
void R_init_tweedie(DllInfo *dll)
{
  // Register C Routines (CEntries), then Fortran Routines (FortranEntries)
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL); 
  
  R_useDynamicSymbols(dll, FALSE); 
}
