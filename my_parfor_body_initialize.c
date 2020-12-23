/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * my_parfor_body_initialize.c
 *
 * Code generation for function 'my_parfor_body_initialize'
 *
 */

/* Include files */
#include "my_parfor_body_initialize.h"
#include "my_parfor_body_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void my_parfor_body_initialize(void)
{
  omp_init_nest_lock(&emlrtNestLockGlobal);
  isInitialized_my_parfor_body_mex = true;
}

/* End of code generation (my_parfor_body_initialize.c) */
