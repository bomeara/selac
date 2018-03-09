#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void initmod_selacHMM();

extern void selacHMM(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"initmod_selacHMM",                   (DL_FUNC) &initmod_selacHMM,                   0},
    {"selacHMM",                           (DL_FUNC) &selacHMM,                           6},
    {NULL, NULL, 0}
};

void R_init_hisse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
