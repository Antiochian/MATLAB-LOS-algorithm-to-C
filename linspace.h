/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.h
 *
 * Code generation for function 'linspace'
 *
 */

#ifndef LINSPACE_H
#define LINSPACE_H

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
  void linspace(double d2, double n, emxArray_real_T *y);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (linspace.h) */
