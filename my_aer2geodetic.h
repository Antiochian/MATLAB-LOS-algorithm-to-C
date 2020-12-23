/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * my_aer2geodetic.h
 *
 * Code generation for function 'my_aer2geodetic'
 *
 */

#ifndef MY_AER2GEODETIC_H
#define MY_AER2GEODETIC_H

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
  void my_aer2geodetic(const emxArray_real_T *az, const emxArray_real_T *elev,
                       double slantRange, double lat0, double lon0, double h0,
                       double spheroid_MeanRadius, double spheroid_Flattening,
                       emxArray_real_T *lat, emxArray_real_T *lon,
                       emxArray_real_T *h);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (my_aer2geodetic.h) */
