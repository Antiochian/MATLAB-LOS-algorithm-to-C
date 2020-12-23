/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mod.h
 *
 * Code generation for function 'mod'
 *
 */

#ifndef MOD_H
#define MOD_H

/* Include files */
#include "my_parfor_body_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void b_mod(const emxArray_real_T *x, emxArray_real_T *r);
  void c_mod(const emxArray_real_T *x, double y, emxArray_real_T *r);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (mod.h) */
