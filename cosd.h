/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cosd.h
 *
 * Code generation for function 'cosd'
 *
 */

#ifndef COSD_H
#define COSD_H

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
  void b_cosd(emxArray_real_T *x);
  void c_cosd(double *x);
  void d_cosd(emxArray_real32_T *x);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (cosd.h) */
