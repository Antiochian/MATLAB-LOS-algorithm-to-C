/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * interp2.h
 *
 * Code generation for function 'interp2'
 *
 */

#ifndef INTERP2_H
#define INTERP2_H

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
  void interp2_local(const emxArray_real_T *V, const emxArray_real_T *Xq, const
                     emxArray_real_T *Yq, emxArray_real_T *Vq);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (interp2.h) */
