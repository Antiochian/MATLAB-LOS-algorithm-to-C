/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.c
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include "linspace.h"
#include "my_parfor_body_emxutil.h"
#include "my_parfor_body_types.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void linspace(double d2, double n, emxArray_real_T *y)
{
  double delta1;
  int k;
  int y_tmp;
  if (!(n >= 0.0)) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    delta1 = floor(n);
    y_tmp = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)delta1;
    emxEnsureCapacity_real_T(y, y_tmp);
    if ((int)delta1 >= 1) {
      y_tmp = (int)delta1 - 1;
      y->data[(int)delta1 - 1] = d2;
      if (y->size[1] >= 2) {
        y->data[0] = 0.0;
        if (y->size[1] >= 3) {
          if ((0.0 == -d2) && ((int)delta1 > 2)) {
            for (k = 2; k <= y_tmp; k++) {
              y->data[k - 1] = d2 * ((double)(((k << 1) - (int)delta1) - 1) /
                ((double)(int)delta1 - 1.0));
            }

            if (((int)delta1 & 1) == 1) {
              y->data[(int)delta1 >> 1] = 0.0;
            }
          } else if ((d2 < 0.0) && (fabs(d2) > 8.9884656743115785E+307)) {
            delta1 = d2 / ((double)y->size[1] - 1.0);
            y_tmp = y->size[1];
            for (k = 0; k <= y_tmp - 3; k++) {
              y->data[k + 1] = delta1 * ((double)k + 1.0);
            }
          } else {
            delta1 = d2 / ((double)y->size[1] - 1.0);
            y_tmp = y->size[1];
            for (k = 0; k <= y_tmp - 3; k++) {
              y->data[k + 1] = ((double)k + 1.0) * delta1;
            }
          }
        }
      }
    }
  }
}

/* End of code generation (linspace.c) */
