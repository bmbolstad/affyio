/*****************************************************
 **
 ** file: init_package.c
 **
 ** Copyright (C) 2013    B. M. Bolstad
 **
 ** aim: Register c code routines so that they can be called in other packages.
 **"
 ** History
 ** May 20, 2013 - Initial version
 **
 *****************************************************/

#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "read_abatch.h"

#if _MSC_VER >= 1000
__declspec(dllexport)
#endif



static const R_CallMethodDef callMethods[]  = {
 {"read_abatch",(DL_FUNC)&read_abatch,7},
  {NULL, NULL, 0}
  };


void R_init_affyio(DllInfo *info){

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);

}
