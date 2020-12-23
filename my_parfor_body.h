/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * my_parfor_body.h
 *
 * Code generation for function 'my_parfor_body'
 *
 */

#ifndef MY_PARFOR_BODY_H
#define MY_PARFOR_BODY_H

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
  extern void my_parfor_body(double t_calculation_steps, double
    t_calculation_step_length, double t_steps, boolean_T use_seasons, const
    emxArray_real_T *decl_arr_local, const emxArray_real_T *h_arr_local, double
    lat, double b_long, double slope, double aspect, double max_distance, double
    height_px, const struct0_T *ref_sphere, const emxArray_real_T
    *extended_elevation_matrix, double extended_raster_ref, double r_moon,
    double sample_dens, double const_decl, double dh, const struct1_T
    *public_raster_struct, emxArray_real_T *theta_matrix);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (my_parfor_body.h) */
