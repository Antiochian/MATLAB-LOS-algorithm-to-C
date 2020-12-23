/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * my_parfor_body.c
 *
 * Code generation for function 'my_parfor_body'
 *
 */

/* Include files */
#include "my_parfor_body.h"
#include "cosd.h"
#include "cotd.h"
#include "mexable_los2.h"
#include "my_aer2geodetic.h"
#include "my_parfor_body_data.h"
#include "my_parfor_body_emxutil.h"
#include "my_parfor_body_initialize.h"
#include "my_parfor_body_types.h"
#include "rt_nonfinite.h"
#include "sind.h"
#include "tand.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void my_parfor_body(double t_calculation_steps, double t_calculation_step_length,
                    double t_steps, boolean_T use_seasons, const emxArray_real_T
                    *decl_arr_local, const emxArray_real_T *h_arr_local, double
                    lat, double b_long, double slope, double aspect, double
                    max_distance, double height_px, const struct0_T *ref_sphere,
                    const emxArray_real_T *extended_elevation_matrix, double
                    extended_raster_ref, double r_moon, double sample_dens,
                    double const_decl, double dh, const struct1_T
                    *public_raster_struct, emxArray_real_T *theta_matrix)
{
  emxArray_boolean_T *sun_vis;
  emxArray_int32_T *r;
  emxArray_int8_T *sigma_ns;
  emxArray_real32_T *sigma_ew;
  emxArray_real32_T *sigma_w;
  emxArray_real_T *az_sun;
  emxArray_real_T *b_z1;
  emxArray_real_T *decl;
  emxArray_real_T *elev_sun;
  emxArray_real_T *h;
  emxArray_real_T *height_sun;
  emxArray_real_T *long_sun;
  emxArray_real_T *theta;
  emxArray_real_T *theta_z;
  emxArray_real_T *z1;
  double apnd;
  double cdiff;
  double ndbl;
  double slant_range;
  double t_arr_start;
  float f;
  int i;
  int k;
  int loop_ub;
  int n;
  int nm1d2;
  int t_calculation_step_idx;
  (void)extended_raster_ref;
  if (!isInitialized_my_parfor_body_mex) {
    my_parfor_body_initialize();
  }

  /*  fprintf('WARNING: Not running MEX file!'); */
  i = theta_matrix->size[0] * theta_matrix->size[1];
  theta_matrix->size[0] = (int)t_calculation_steps;
  theta_matrix->size[1] = (int)t_calculation_step_length;
  emxEnsureCapacity_real_T(theta_matrix, i);
  loop_ub = (int)t_calculation_steps * (int)t_calculation_step_length;
  for (i = 0; i < loop_ub; i++) {
    theta_matrix->data[i] = 0.0;
  }

  loop_ub = (int)t_calculation_steps - 1;

#pragma omp parallel \
 num_threads(omp_get_max_threads()) \
 private(r,z1,b_z1,theta,decl,sigma_ns,sigma_w,theta_z,sigma_ew,az_sun,elev_sun,h,long_sun,height_sun,sun_vis,slant_range,t_arr_start,k,nm1d2,n,f,ndbl,apnd,cdiff)

  {
    emxInit_int32_T(&r, 2);
    emxInit_real_T(&z1, 2);
    emxInit_real_T(&b_z1, 2);
    emxInit_real_T(&theta, 2);
    emxInit_real_T(&decl, 2);
    emxInit_int8_T(&sigma_ns, 2);
    emxInit_real32_T(&sigma_w, 2);
    emxInit_real_T(&theta_z, 2);
    emxInit_real32_T(&sigma_ew, 2);
    emxInit_real_T(&az_sun, 2);
    emxInit_real_T(&elev_sun, 2);
    emxInit_real_T(&h, 2);
    emxInit_real_T(&long_sun, 2);
    emxInit_real_T(&height_sun, 2);
    emxInit_boolean_T(&sun_vis, 2);

#pragma omp for nowait

    for (t_calculation_step_idx = 0; t_calculation_step_idx <= loop_ub;
         t_calculation_step_idx++) {
      /* parfor */
      k = theta->size[0] * theta->size[1];
      theta->size[0] = 1;
      theta->size[1] = (int)t_calculation_step_length;
      emxEnsureCapacity_real_T(theta, k);
      nm1d2 = (int)t_calculation_step_length;
      for (k = 0; k < nm1d2; k++) {
        theta->data[k] = rtNaN;
      }

      t_arr_start = (((double)t_calculation_step_idx + 1.0) - 1.0) *
        t_calculation_step_length + 1.0;
      slant_range = (t_arr_start + t_calculation_step_length) - 1.0;
      if (slant_range > t_steps) {
        slant_range = t_steps;
      }

      if (rtIsNaN(t_arr_start) || rtIsNaN(slant_range)) {
        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        elev_sun->size[1] = 1;
        emxEnsureCapacity_real_T(elev_sun, k);
        elev_sun->data[0] = rtNaN;
      } else if (slant_range < t_arr_start) {
        elev_sun->size[0] = 1;
        elev_sun->size[1] = 0;
      } else if ((rtIsInf(t_arr_start) || rtIsInf(slant_range)) && (t_arr_start ==
                  slant_range)) {
        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        elev_sun->size[1] = 1;
        emxEnsureCapacity_real_T(elev_sun, k);
        elev_sun->data[0] = rtNaN;
      } else if (t_arr_start == t_arr_start) {
        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        nm1d2 = (int)floor(slant_range - t_arr_start);
        elev_sun->size[1] = nm1d2 + 1;
        emxEnsureCapacity_real_T(elev_sun, k);
        for (k = 0; k <= nm1d2; k++) {
          elev_sun->data[k] = t_arr_start + (double)k;
        }
      } else {
        ndbl = floor((slant_range - t_arr_start) + 0.5);
        apnd = t_arr_start + ndbl;
        cdiff = apnd - slant_range;
        if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(t_arr_start, fabs
             (slant_range))) {
          ndbl++;
          apnd = slant_range;
        } else if (cdiff > 0.0) {
          apnd = t_arr_start + (ndbl - 1.0);
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        elev_sun->size[1] = n;
        emxEnsureCapacity_real_T(elev_sun, k);
        if (n > 0) {
          elev_sun->data[0] = t_arr_start;
          if (n > 1) {
            elev_sun->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 0; k <= nm1d2 - 2; k++) {
              elev_sun->data[k + 1] = t_arr_start + ((double)k + 1.0);
              elev_sun->data[(n - k) - 2] = apnd - ((double)k + 1.0);
            }

            if (nm1d2 << 1 == n - 1) {
              elev_sun->data[nm1d2] = (t_arr_start + apnd) / 2.0;
            } else {
              elev_sun->data[nm1d2] = t_arr_start + (double)nm1d2;
              elev_sun->data[nm1d2 + 1] = apnd - (double)nm1d2;
            }
          }
        }
      }

      /*  array of time values to calculate for at this iteration */
      if (elev_sun->size[1] != 0) {
        if (use_seasons) {
          k = decl->size[0] * decl->size[1];
          decl->size[0] = 1;
          decl->size[1] = elev_sun->size[1];
          emxEnsureCapacity_real_T(decl, k);
          nm1d2 = elev_sun->size[0] * elev_sun->size[1];
          for (k = 0; k < nm1d2; k++) {
            decl->data[k] = decl_arr_local->data[(int)elev_sun->data[k] - 1];
          }

          k = h->size[0] * h->size[1];
          h->size[0] = 1;
          h->size[1] = elev_sun->size[1];
          emxEnsureCapacity_real_T(h, k);
          nm1d2 = elev_sun->size[0] * elev_sun->size[1];
          for (k = 0; k < nm1d2; k++) {
            h->data[k] = h_arr_local->data[(int)elev_sun->data[k] - 1];
          }
        } else {
          k = decl->size[0] * decl->size[1];
          decl->size[0] = 1;
          decl->size[1] = elev_sun->size[1];
          emxEnsureCapacity_real_T(decl, k);
          nm1d2 = elev_sun->size[1];
          for (k = 0; k < nm1d2; k++) {
            decl->data[k] = const_decl;
          }

          k = h->size[0] * h->size[1];
          h->size[0] = 1;
          h->size[1] = elev_sun->size[1];
          emxEnsureCapacity_real_T(h, k);
          nm1d2 = elev_sun->size[0] * elev_sun->size[1];
          for (k = 0; k < nm1d2; k++) {
            h->data[k] = dh * elev_sun->data[k] - dh;
          }
        }

        k = h->size[0] * h->size[1];
        n = h->size[0] * h->size[1];
        h->size[0] = 1;
        emxEnsureCapacity_real_T(h, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          h->data[k] += 180.0;
        }

        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = h->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = h->size[1];
        for (k = 0; k < nm1d2; k++) {
          t_arr_start = h->data[k];
          if (rtIsNaN(t_arr_start) || rtIsInf(t_arr_start)) {
            slant_range = rtNaN;
          } else if (t_arr_start == 0.0) {
            slant_range = 0.0;
          } else {
            slant_range = fmod(t_arr_start, 360.0);
            if (slant_range == 0.0) {
              slant_range = 0.0;
            } else {
              if (t_arr_start < 0.0) {
                slant_range += 360.0;
              }
            }
          }

          height_sun->data[k] = slant_range;
        }

        k = h->size[0] * h->size[1];
        h->size[0] = 1;
        h->size[1] = height_sun->size[1];
        emxEnsureCapacity_real_T(h, k);
        nm1d2 = height_sun->size[0] * height_sun->size[1];
        for (k = 0; k < nm1d2; k++) {
          h->data[k] = height_sun->data[k] - 180.0;
        }

        k = sigma_ew->size[0] * sigma_ew->size[1];
        sigma_ew->size[0] = 1;
        sigma_ew->size[1] = elev_sun->size[1];
        emxEnsureCapacity_real32_T(sigma_ew, k);
        nm1d2 = elev_sun->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_ew->data[k] = 1.0F;
        }

        nm1d2 = h->size[1];
        k = long_sun->size[0] * long_sun->size[1];
        long_sun->size[0] = 1;
        long_sun->size[1] = h->size[1];
        emxEnsureCapacity_real_T(long_sun, k);
        for (k = 0; k < nm1d2; k++) {
          long_sun->data[k] = fabs(h->data[k]);
        }

        t_arr_start = lat;
        b_cotd(&t_arr_start);
        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = decl->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = decl->size[0] * decl->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = decl->data[k];
        }

        b_tand(height_sun);
        k = height_sun->size[0] * height_sun->size[1];
        n = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        emxEnsureCapacity_real_T(height_sun, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          height_sun->data[k] *= t_arr_start;
        }

        k = b_z1->size[0] * b_z1->size[1];
        b_z1->size[0] = 1;
        n = height_sun->size[1];
        b_z1->size[1] = height_sun->size[1];
        emxEnsureCapacity_real_T(b_z1, k);
        for (k = 0; k < n; k++) {
          b_z1->data[k] = fmax(height_sun->data[k], -1.0);
        }

        k = z1->size[0] * z1->size[1];
        z1->size[0] = 1;
        n = b_z1->size[1];
        z1->size[1] = b_z1->size[1];
        emxEnsureCapacity_real_T(z1, k);
        for (k = 0; k < n; k++) {
          z1->data[k] = fmin(b_z1->data[k], 1.0);
        }

        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = z1->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = z1->size[0] * z1->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = z1->data[k];
        }

        nm1d2 = z1->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = 57.295779513082323 * acos(height_sun->data[k]);
        }

        n = long_sun->size[1];
        for (k = 0; k < n; k++) {
          if (long_sun->data[k] > height_sun->data[k]) {
            sigma_ew->data[k] = -1.0F;
          }
        }

        k = sigma_ns->size[0] * sigma_ns->size[1];
        sigma_ns->size[0] = 1;
        sigma_ns->size[1] = elev_sun->size[1];
        emxEnsureCapacity_int8_T(sigma_ns, k);
        nm1d2 = elev_sun->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_ns->data[k] = 1;
        }

        k = sun_vis->size[0] * sun_vis->size[1];
        sun_vis->size[0] = 1;
        sun_vis->size[1] = decl->size[1];
        emxEnsureCapacity_boolean_T(sun_vis, k);
        nm1d2 = decl->size[0] * decl->size[1];
        for (k = 0; k < nm1d2; k++) {
          sun_vis->data[k] = (lat * (lat - decl->data[k]) < 0.0);
        }

        n = sun_vis->size[1] - 1;
        nm1d2 = 0;
        for (k = 0; k <= n; k++) {
          if (sun_vis->data[k]) {
            nm1d2++;
          }
        }

        k = r->size[0] * r->size[1];
        r->size[0] = 1;
        r->size[1] = nm1d2;
        emxEnsureCapacity_int32_T(r, k);
        nm1d2 = 0;
        for (k = 0; k <= n; k++) {
          if (sun_vis->data[k]) {
            r->data[nm1d2] = k + 1;
            nm1d2++;
          }
        }

        nm1d2 = r->size[0] * r->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_ns->data[r->data[k] - 1] = -1;
        }

        k = sigma_w->size[0] * sigma_w->size[1];
        sigma_w->size[0] = 1;
        sigma_w->size[1] = elev_sun->size[1];
        emxEnsureCapacity_real32_T(sigma_w, k);
        nm1d2 = elev_sun->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_w->data[k] = 1.0F;
        }

        n = h->size[1];
        for (k = 0; k < n; k++) {
          if (h->data[k] < 0.0) {
            sigma_w->data[k] = -1.0F;
          }
        }

        k = long_sun->size[0] * long_sun->size[1];
        long_sun->size[0] = 1;
        long_sun->size[1] = decl->size[1];
        emxEnsureCapacity_real_T(long_sun, k);
        nm1d2 = decl->size[0] * decl->size[1];
        for (k = 0; k < nm1d2; k++) {
          long_sun->data[k] = decl->data[k];
        }

        b_cosd(long_sun);
        k = theta_z->size[0] * theta_z->size[1];
        theta_z->size[0] = 1;
        theta_z->size[1] = decl->size[1];
        emxEnsureCapacity_real_T(theta_z, k);
        nm1d2 = decl->size[0] * decl->size[1];
        for (k = 0; k < nm1d2; k++) {
          theta_z->data[k] = decl->data[k];
        }

        b_sind(theta_z);
        t_arr_start = lat;
        c_sind(&t_arr_start);
        slant_range = lat;
        c_cosd(&slant_range);
        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = h->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = h->size[0] * h->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = h->data[k];
        }

        b_cosd(height_sun);
        k = theta_z->size[0] * theta_z->size[1];
        n = theta_z->size[0] * theta_z->size[1];
        theta_z->size[0] = 1;
        emxEnsureCapacity_real_T(theta_z, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          theta_z->data[k] = theta_z->data[k] * t_arr_start + long_sun->data[k] *
            slant_range * height_sun->data[k];
        }

        nm1d2 = theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          theta_z->data[k] = acos(theta_z->data[k]);
        }

        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = theta_z->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = theta_z->size[0] * theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = theta_z->data[k];
        }

        nm1d2 = theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = sin(height_sun->data[k]);
        }

        b_sind(h);
        k = h->size[0] * h->size[1];
        n = h->size[0] * h->size[1];
        h->size[0] = 1;
        emxEnsureCapacity_real_T(h, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          h->data[k] = h->data[k] * long_sun->data[k] / height_sun->data[k];
        }

        nm1d2 = h->size[1];
        for (k = 0; k < nm1d2; k++) {
          h->data[k] = 57.295779513082323 * asin(h->data[k]);
        }

        n = theta_z->size[1];
        for (k = 0; k < n; k++) {
          if (theta_z->data[k] == 0.0) {
            h->data[k] = 0.0;
          }
        }

        /*  Avoid division by 0 error */
        /*  gamma_so may end up complex if it is ~90deg as rounding errors in */
        /*  the asind(...) function can make its argument slightly larger */
        /*  than 1, making gamma_so slightly complex. Therefore, only use the */
        /*  real part of gamma_so to avoid the tiny complex part affecting */
        /*  later calculations. */
        k = sigma_ew->size[0] * sigma_ew->size[1];
        n = sigma_ew->size[0] * sigma_ew->size[1];
        sigma_ew->size[0] = 1;
        emxEnsureCapacity_real32_T(sigma_ew, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          sigma_ew->data[k] *= (float)sigma_ns->data[k];
        }

        k = sigma_ew->size[0] * sigma_ew->size[1];
        n = sigma_ew->size[0] * sigma_ew->size[1];
        sigma_ew->size[0] = 1;
        emxEnsureCapacity_real32_T(sigma_ew, n);
        nm1d2 = k - 1;
        for (k = 0; k <= nm1d2; k++) {
          f = sigma_ew->data[k];
          f = f * (float)h->data[k] + (1.0F - f) / 2.0F * sigma_w->data[k] *
            180.0F;
          sigma_ew->data[k] = f;
        }

        k = height_sun->size[0] * height_sun->size[1];
        height_sun->size[0] = 1;
        height_sun->size[1] = theta_z->size[1];
        emxEnsureCapacity_real_T(height_sun, k);
        nm1d2 = theta_z->size[0] * theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = theta_z->data[k];
        }

        nm1d2 = theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          height_sun->data[k] = sin(height_sun->data[k]);
        }

        k = long_sun->size[0] * long_sun->size[1];
        long_sun->size[0] = 1;
        long_sun->size[1] = theta_z->size[1];
        emxEnsureCapacity_real_T(long_sun, k);
        nm1d2 = theta_z->size[0] * theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          long_sun->data[k] = theta_z->data[k];
        }

        nm1d2 = theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          long_sun->data[k] = cos(long_sun->data[k]);
        }

        t_arr_start = slope;
        c_cosd(&t_arr_start);
        slant_range = slope;
        c_sind(&slant_range);
        k = sigma_w->size[0] * sigma_w->size[1];
        sigma_w->size[0] = 1;
        sigma_w->size[1] = sigma_ew->size[1];
        emxEnsureCapacity_real32_T(sigma_w, k);
        nm1d2 = sigma_ew->size[0] * sigma_ew->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_w->data[k] = sigma_ew->data[k] - (float)aspect;
        }

        d_cosd(sigma_w);
        k = sigma_w->size[0] * sigma_w->size[1];
        sigma_w->size[0] = 1;
        sigma_w->size[1] = long_sun->size[1];
        emxEnsureCapacity_real32_T(sigma_w, k);
        nm1d2 = long_sun->size[0] * long_sun->size[1] - 1;
        for (k = 0; k <= nm1d2; k++) {
          sigma_w->data[k] = (float)(long_sun->data[k] * t_arr_start) + (float)
            (height_sun->data[k] * slant_range) * sigma_w->data[k];
        }

        nm1d2 = sigma_w->size[1];
        for (k = 0; k < nm1d2; k++) {
          sigma_w->data[k] = acosf(sigma_w->data[k]);
        }

        nm1d2 = sigma_w->size[1];
        for (k = 0; k < nm1d2; k++) {
          theta->data[k] = sigma_w->data[k];
        }

        n = theta_z->size[1];
        for (k = 0; k < n; k++) {
          if (theta_z->data[k] > 1.5707963267948966) {
            theta->data[k] = rtNaN;
          }
        }

        n = theta->size[1];
        for (k = 0; k < n; k++) {
          if (theta->data[k] > 1.5707963267948966) {
            theta->data[k] = rtNaN;
          }
        }

        k = az_sun->size[0] * az_sun->size[1];
        az_sun->size[0] = 1;
        az_sun->size[1] = elev_sun->size[1];
        emxEnsureCapacity_real_T(az_sun, k);
        nm1d2 = elev_sun->size[1];
        for (k = 0; k < nm1d2; k++) {
          az_sun->data[k] = rtNaN;
        }

        n = decl->size[1];
        for (k = 0; k < n; k++) {
          if (lat >= decl->data[k]) {
            az_sun->data[k] = sigma_ew->data[k] + 180.0F;
          }
        }

        n = decl->size[1];
        for (k = 0; k < n; k++) {
          if (lat < decl->data[k]) {
            az_sun->data[k] = -sigma_ew->data[k];
          }
        }

        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        elev_sun->size[1] = theta_z->size[1];
        emxEnsureCapacity_real_T(elev_sun, k);
        nm1d2 = theta_z->size[0] * theta_z->size[1];
        for (k = 0; k < nm1d2; k++) {
          elev_sun->data[k] = 90.0 - 57.295779513082323 * theta_z->data[k];
        }

        slant_range = 1.1 * max_distance;

        /*  ensure outside of grid */
        /*  switch back up to double precision here so that built in */
        /*  matlab functions work */
        my_aer2geodetic(az_sun, elev_sun, slant_range, lat, b_long, height_px,
                        ref_sphere->MeanRadius, ref_sphere->Flattening, h,
                        long_sun, height_sun);

        /*          sun_vis1 = my_los2(extended_elevation_matrix, extended_raster_ref, lat*ones(size(lat_sun)), long*ones(size(lat_sun)), lat_sun, long_sun, 0, height_sun, 'AGL', 'MSL', r_moon); */
        k = elev_sun->size[0] * elev_sun->size[1];
        elev_sun->size[0] = 1;
        elev_sun->size[1] = h->size[1];
        emxEnsureCapacity_real_T(elev_sun, k);
        nm1d2 = h->size[1];
        for (k = 0; k < nm1d2; k++) {
          elev_sun->data[k] = lat;
        }

        k = az_sun->size[0] * az_sun->size[1];
        az_sun->size[0] = 1;
        az_sun->size[1] = h->size[1];
        emxEnsureCapacity_real_T(az_sun, k);
        nm1d2 = h->size[1];
        for (k = 0; k < nm1d2; k++) {
          az_sun->data[k] = b_long;
        }

        mexable_los2(extended_elevation_matrix, elev_sun, az_sun, h, long_sun,
                     height_sun, r_moon, sample_dens, public_raster_struct,
                     sun_vis);
        n = theta_z->size[1];
        for (k = 0; k < n; k++) {
          if (theta_z->data[k] == 0.0) {
            sun_vis->data[k] = true;
          }
        }

        n = sun_vis->size[1];
        for (k = 0; k < n; k++) {
          if (!sun_vis->data[k]) {
            theta->data[k] = rtNaN;
          }
        }

        if (theta->size[1] < t_calculation_step_length) {
          if ((double)theta->size[1] + 1.0 > t_calculation_step_length) {
            k = 0;
            n = 0;
          } else {
            k = theta->size[1];
            n = (int)t_calculation_step_length;
          }

          nm1d2 = n - k;
          for (n = 0; n < nm1d2; n++) {
            theta->data[k + n] = rtNaN;
          }
        }

        nm1d2 = theta->size[1];
        for (k = 0; k < nm1d2; k++) {
          theta_matrix->data[t_calculation_step_idx + theta_matrix->size[0] * k]
            = theta->data[k];
        }
      }
    }

    emxFree_boolean_T(&sun_vis);
    emxFree_real_T(&height_sun);
    emxFree_real_T(&long_sun);
    emxFree_real_T(&h);
    emxFree_real_T(&elev_sun);
    emxFree_real_T(&az_sun);
    emxFree_real32_T(&sigma_ew);
    emxFree_real_T(&theta_z);
    emxFree_real32_T(&sigma_w);
    emxFree_int8_T(&sigma_ns);
    emxFree_real_T(&decl);
    emxFree_real_T(&theta);
    emxFree_real_T(&b_z1);
    emxFree_real_T(&z1);
    emxFree_int32_T(&r);
  }
}

/* End of code generation (my_parfor_body.c) */
