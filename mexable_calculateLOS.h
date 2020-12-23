/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mexable_calculateLOS.h
 *
 * Code generation for function 'mexable_calculateLOS'
 *
 */

#ifndef MEXABLE_CALCULATELOS_H
#define MEXABLE_CALCULATELOS_H

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
  void mexable_calculateLOS(const emxArray_real_T *Z, double lat1, double lon1,
    double lat2, double lon2, double talt, double actualradius, double
    apparentradius, double sample_dens, const struct1_T *public_raster_struct,
    emxArray_boolean_T *vis0, emxArray_real_T *dist, emxArray_real_T *h,
    emxArray_real_T *lattrk, emxArray_real_T *lontrk, emxArray_real_T *x1,
    emxArray_real_T *z1, emxArray_real_T *x2, emxArray_real_T *z2);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (mexable_calculateLOS.h) */
