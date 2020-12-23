/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cotd.c
 *
 * Code generation for function 'cotd'
 *
 */

/* Include files */
#include "cotd.h"
#include "my_parfor_body_rtwutil.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void b_cotd(double *x)
{
  double absx;
  double z;
  signed char n;
  if (rtIsInf(*x) || rtIsNaN(*x)) {
    absx = rtNaN;
  } else {
    *x = rt_remd_snf(*x, 360.0);
    absx = fabs(*x);
    if (absx > 180.0) {
      if (*x > 0.0) {
        *x -= 360.0;
      } else {
        *x += 360.0;
      }

      absx = fabs(*x);
    }

    if (absx <= 45.0) {
      *x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (*x > 0.0) {
        *x = 0.017453292519943295 * (*x - 90.0);
        n = 1;
      } else {
        *x = 0.017453292519943295 * (*x + 90.0);
        n = -1;
      }
    } else if (*x > 0.0) {
      *x = 0.017453292519943295 * (*x - 180.0);
      n = 2;
    } else {
      *x = 0.017453292519943295 * (*x + 180.0);
      n = -2;
    }

    absx = tan(*x);
    if ((n == 1) || (n == -1)) {
      z = 1.0 / absx;
      absx = -(1.0 / absx);
      if (rtIsInf(absx) && (n == 1)) {
        absx = z;
      }
    }
  }

  *x = 1.0 / absx;
}

/* End of code generation (cotd.c) */
