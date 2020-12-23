/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mod.c
 *
 * Code generation for function 'mod'
 *
 */

/* Include files */
#include "mod.h"
#include "my_parfor_body_emxutil.h"
#include "my_parfor_body_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void b_mod(const emxArray_real_T *x, emxArray_real_T *r)
{
  double b_r;
  double q;
  int k;
  int nx;
  boolean_T rEQ0;
  nx = r->size[0];
  r->size[0] = x->size[0];
  emxEnsureCapacity_real_T(r, nx);
  nx = x->size[0];
  for (k = 0; k < nx; k++) {
    if (rtIsNaN(x->data[k]) || rtIsInf(x->data[k])) {
      b_r = rtNaN;
    } else if (x->data[k] == 0.0) {
      b_r = 0.0;
    } else {
      b_r = fmod(x->data[k], 6.2831853071795862);
      rEQ0 = (b_r == 0.0);
      if (!rEQ0) {
        q = fabs(x->data[k] / 6.2831853071795862);
        rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
      }

      if (rEQ0) {
        b_r = 0.0;
      } else {
        if (x->data[k] < 0.0) {
          b_r += 6.2831853071795862;
        }
      }
    }

    r->data[k] = b_r;
  }
}

void c_mod(const emxArray_real_T *x, double y, emxArray_real_T *r)
{
  double b_r;
  double q;
  int k;
  int nx;
  boolean_T rEQ0;
  nx = r->size[0];
  r->size[0] = x->size[0];
  emxEnsureCapacity_real_T(r, nx);
  nx = x->size[0];
  for (k = 0; k < nx; k++) {
    b_r = x->data[k];
    if (y == 0.0) {
      if (x->data[k] == 0.0) {
        b_r = y;
      }
    } else if (rtIsNaN(x->data[k]) || rtIsNaN(y) || rtIsInf(x->data[k])) {
      b_r = rtNaN;
    } else if (x->data[k] == 0.0) {
      b_r = 0.0 / y;
    } else if (rtIsInf(y)) {
      if ((y < 0.0) != (x->data[k] < 0.0)) {
        b_r = y;
      }
    } else {
      b_r = fmod(x->data[k], y);
      rEQ0 = (b_r == 0.0);
      if ((!rEQ0) && (y > floor(y))) {
        q = fabs(x->data[k] / y);
        rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
      }

      if (rEQ0) {
        b_r = y * 0.0;
      } else {
        if ((x->data[k] < 0.0) != (y < 0.0)) {
          b_r += y;
        }
      }
    }

    r->data[k] = b_r;
  }
}

/* End of code generation (mod.c) */
