#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP kzp2w(SEXP, SEXP, SEXP);
extern SEXP kzpg(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"kzp2w", (DL_FUNC) &kzp2w, 3},
    {"kzpg",  (DL_FUNC) &kzpg,  3},
    {NULL, NULL, 0}
};

void R_init_kzfs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
