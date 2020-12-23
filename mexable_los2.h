/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mexable_los2.h
 *
 * Code generation for function 'mexable_los2'
 *
 */

#ifndef MEXABLE_LOS2_H
#define MEXABLE_LOS2_H

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
  void mexable_los2(const emxArray_real_T *Z, const emxArray_real_T *lat1, const
                    emxArray_real_T *lon1, const emxArray_real_T *lat2, const
                    emxArray_real_T *lon2, const emxArray_real_T *inp_talt,
                    double actualradius, double sample_dens, const struct1_T
                    *public_raster_struct, emxArray_boolean_T *vis);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (mexable_los2.h) */
