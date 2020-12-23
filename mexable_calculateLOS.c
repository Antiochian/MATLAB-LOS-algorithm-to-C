/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mexable_calculateLOS.c
 *
 * Code generation for function 'mexable_calculateLOS'
 *
 */

/* Include files */
#include "mexable_calculateLOS.h"
#include "interp2.h"
#include "linspace.h"
#include "mod.h"
#include "my_parfor_body_emxutil.h"
#include "my_parfor_body_rtwutil.h"
#include "my_parfor_body_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Declarations */
static void adjustterrain(const emxArray_real_T *arclen, const emxArray_real_T
  *zin, double apparentradius, emxArray_real_T *x, emxArray_real_T *z);

/* Function Definitions */
static void adjustterrain(const emxArray_real_T *arclen, const emxArray_real_T
  *zin, double apparentradius, emxArray_real_T *x, emxArray_real_T *z)
{
  emxArray_real_T *phi;
  emxArray_real_T *phi0;
  int i;
  int idx;
  int j;

  /* ----------------------------------------------------------------------- */
  /*  Adjust the terrain slice for the curvature of the sphere. The radius */
  /*  may potentially be different from the actual body, for example to */
  /*  model refraction of radio waves. */
  i = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = zin->size[1];
  emxEnsureCapacity_real_T(z, i);
  idx = zin->size[0] * zin->size[1];
  for (i = 0; i < idx; i++) {
    z->data[i] = apparentradius + zin->data[i];
  }

  emxInit_real_T(&phi0, 2);
  i = phi0->size[0] * phi0->size[1];
  phi0->size[0] = 1;
  phi0->size[1] = arclen->size[1];
  emxEnsureCapacity_real_T(phi0, i);
  idx = arclen->size[0] * arclen->size[1];
  for (i = 0; i < idx; i++) {
    phi0->data[i] = arclen->data[i] / apparentradius;
  }

  emxInit_real_T(&phi, 2);
  if (phi0->size[1] == 1) {
    i = phi->size[0] * phi->size[1];
    phi->size[0] = 1;
    phi->size[1] = z->size[1];
    emxEnsureCapacity_real_T(phi, i);
    idx = -1;
    i = z->size[1];
    for (j = 0; j < i; j++) {
      idx++;
      phi->data[idx] = phi0->data[0];
    }
  } else {
    i = phi->size[0] * phi->size[1];
    phi->size[0] = 1;
    phi->size[1] = phi0->size[1];
    emxEnsureCapacity_real_T(phi, i);
    idx = phi0->size[0] * phi0->size[1];
    for (i = 0; i < idx; i++) {
      phi->data[i] = phi0->data[i];
    }
  }

  emxFree_real_T(&phi0);
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = phi->size[1];
  emxEnsureCapacity_real_T(x, i);
  idx = phi->size[0] * phi->size[1];
  for (i = 0; i < idx; i++) {
    x->data[i] = phi->data[i];
  }

  idx = phi->size[1];
  for (j = 0; j < idx; j++) {
    x->data[j] = sin(x->data[j]);
  }

  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = z->size[1];
  emxEnsureCapacity_real_T(x, i);
  idx = z->size[0] * z->size[1] - 1;
  for (i = 0; i <= idx; i++) {
    x->data[i] *= z->data[i];
  }

  idx = phi->size[1];
  for (j = 0; j < idx; j++) {
    phi->data[j] = cos(phi->data[j]);
  }

  i = z->size[0] * z->size[1];
  idx = z->size[0] * z->size[1];
  z->size[0] = 1;
  emxEnsureCapacity_real_T(z, idx);
  idx = i - 1;
  for (i = 0; i <= idx; i++) {
    z->data[i] = z->data[i] * phi->data[i] - apparentradius;
  }

  emxFree_real_T(&phi);
}

void mexable_calculateLOS(const emxArray_real_T *Z, double lat1, double lon1,
  double lat2, double lon2, double talt, double actualradius, double
  apparentradius, double sample_dens, const struct1_T *public_raster_struct,
  emxArray_boolean_T *vis0, emxArray_real_T *dist, emxArray_real_T *h,
  emxArray_real_T *lattrk, emxArray_real_T *lontrk, emxArray_real_T *x1,
  emxArray_real_T *z1, emxArray_real_T *x2, emxArray_real_T *z2)
{
  emxArray_boolean_T *positiveInput;
  emxArray_boolean_T *r;
  emxArray_boolean_T *r1;
  emxArray_boolean_T *r4;
  emxArray_boolean_T *vis;
  emxArray_int32_T *r2;
  emxArray_int32_T *r3;
  emxArray_int32_T *r5;
  emxArray_real_T *ang;
  emxArray_real_T *ang2;
  emxArray_real_T *b_dist;
  emxArray_real_T *cosDelta;
  emxArray_real_T *sinDelta;
  emxArray_real_T *x;
  emxArray_real_T *xi;
  emxArray_real_T *y;
  double a;
  double a_tmp;
  double az_tmp;
  double b_a;
  double cosAz;
  double cosphi1_tmp;
  double cosphi2;
  double lambda1;
  double maxdist;
  double spacingInDegrees;
  double w;
  int i;
  int k;
  int nx;
  boolean_T b;

  /*  Perform the line-of-sight computations needed by LOS2 and VIEWSHED. */
  /*  Copyright 2014-2015 The MathWorks, Inc. */
  /*  Sample at slightly less than the elevation grid spacing. */
  spacingInDegrees = 0.9 / sample_dens;

  /*  Intermediate points along the great circle arc between start and end. */
  /* -------------------------------------------------------------------------- */
  /*   Compute sort of maximum angular distance between the end points. */
  maxdist = fmax(fabs(lat2 - lat1), fabs(lon2 - lon1));
  emxInit_real_T(&xi, 1);
  emxInit_real_T(&ang, 2);
  emxInit_boolean_T(&r, 1);
  emxInit_boolean_T(&r1, 1);
  emxInit_real_T(&cosDelta, 1);
  emxInit_real_T(&sinDelta, 1);
  emxInit_real_T(&y, 1);
  if (maxdist > spacingInDegrees) {
    /*   Insert points using linear interpolation. */
    w = 0.017453292519943295 * lat1;
    lambda1 = 0.017453292519943295 * lon1;
    cosAz = 0.017453292519943295 * lat2;

    /* -------------------------------------------------------------------------- */
    /*  Interpolate regularly spaced points along a great circle. */
    /* -------------------------------------------------------------------------- */
    /*  Great circle distance and azimuth between points on a sphere, using the */
    /*  Haversine Formula for distance.  All angles are in radians. */
    cosphi1_tmp = cos(w);
    cosphi2 = cos(cosAz);
    a = sin((cosAz - w) / 2.0);
    a_tmp = 0.017453292519943295 * lon2 - lambda1;
    b_a = sin(a_tmp / 2.0);
    az_tmp = sin(w);
    w = rt_atan2d_snf(cosphi2 * sin(a_tmp), cosphi1_tmp * sin(cosAz) - az_tmp *
                      cosphi2 * cos(a_tmp));
    linspace(2.0 * asin(sqrt(a * a + cosphi1_tmp * cosphi2 * (b_a * b_a))), ceil
             (maxdist / spacingInDegrees) + 1.0, ang);
    i = dist->size[0];
    dist->size[0] = ang->size[1];
    emxEnsureCapacity_real_T(dist, i);
    nx = ang->size[1];
    for (i = 0; i < nx; i++) {
      dist->data[i] = ang->data[i];
    }

    /* -------------------------------------------------------------------------- */
    /*  Points on a great circles given specified start point, azimuths and */
    /*  spherical distances.  All angles are in radians. */
    /*  Reference */
    /*  --------- */
    /*  J. P. Snyder, "Map Projections - A Working Manual,"  US Geological Survey */
    /*  Professional Paper 1395, US Government Printing Office, Washington, DC, */
    /*  1987, pp. 29-32. */
    cosAz = cos(w);
    i = cosDelta->size[0];
    cosDelta->size[0] = dist->size[0];
    emxEnsureCapacity_real_T(cosDelta, i);
    nx = dist->size[0];
    for (i = 0; i < nx; i++) {
      cosDelta->data[i] = dist->data[i];
    }

    nx = dist->size[0];
    for (k = 0; k < nx; k++) {
      cosDelta->data[k] = cos(cosDelta->data[k]);
    }

    i = sinDelta->size[0];
    sinDelta->size[0] = dist->size[0];
    emxEnsureCapacity_real_T(sinDelta, i);
    nx = dist->size[0];
    for (i = 0; i < nx; i++) {
      sinDelta->data[i] = dist->data[i];
    }

    nx = dist->size[0];
    for (k = 0; k < nx; k++) {
      sinDelta->data[k] = sin(sinDelta->data[k]);
    }

    w = sin(w);
    i = y->size[0];
    y->size[0] = sinDelta->size[0];
    emxEnsureCapacity_real_T(y, i);
    nx = sinDelta->size[0];
    for (i = 0; i < nx; i++) {
      y->data[i] = sinDelta->data[i] * w;
    }

    emxInit_real_T(&x, 1);
    i = x->size[0];
    x->size[0] = cosDelta->size[0];
    emxEnsureCapacity_real_T(x, i);
    nx = cosDelta->size[0];
    for (i = 0; i < nx; i++) {
      x->data[i] = cosphi1_tmp * cosDelta->data[i] - az_tmp * sinDelta->data[i] *
        cosAz;
    }

    i = xi->size[0];
    if (y->size[0] <= x->size[0]) {
      xi->size[0] = y->size[0];
    } else {
      xi->size[0] = x->size[0];
    }

    emxEnsureCapacity_real_T(xi, i);
    if (y->size[0] <= x->size[0]) {
      nx = y->size[0];
    } else {
      nx = x->size[0];
    }

    for (k = 0; k < nx; k++) {
      xi->data[k] = rt_atan2d_snf(y->data[k], x->data[k]);
    }

    nx = xi->size[0];
    for (i = 0; i < nx; i++) {
      xi->data[i] += lambda1;
    }

    /*  lambdaTrack = wrapToPi(lambdaTrack); */
    i = r->size[0];
    r->size[0] = xi->size[0];
    emxEnsureCapacity_boolean_T(r, i);
    nx = xi->size[0];
    for (i = 0; i < nx; i++) {
      r->data[i] = (xi->data[i] < -3.1415926535897931);
    }

    i = r1->size[0];
    r1->size[0] = xi->size[0];
    emxEnsureCapacity_boolean_T(r1, i);
    nx = xi->size[0];
    for (i = 0; i < nx; i++) {
      r1->data[i] = (3.1415926535897931 < xi->data[i]);
    }

    k = r->size[0] - 1;
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (r->data[i] || r1->data[i]) {
        nx++;
      }
    }

    emxInit_int32_T(&r2, 1);
    i = r2->size[0];
    r2->size[0] = nx;
    emxEnsureCapacity_int32_T(r2, i);
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (r->data[i] || r1->data[i]) {
        r2->data[nx] = i + 1;
        nx++;
      }
    }

    i = y->size[0];
    y->size[0] = r2->size[0];
    emxEnsureCapacity_real_T(y, i);
    nx = r2->size[0];
    for (i = 0; i < nx; i++) {
      y->data[i] = xi->data[r2->data[i] - 1] + 3.1415926535897931;
    }

    emxFree_int32_T(&r2);
    emxInit_boolean_T(&positiveInput, 1);

    /* -------------------------------------------------------------------------- */
    i = positiveInput->size[0];
    positiveInput->size[0] = y->size[0];
    emxEnsureCapacity_boolean_T(positiveInput, i);
    nx = y->size[0];
    for (i = 0; i < nx; i++) {
      positiveInput->data[i] = (y->data[i] > 0.0);
    }

    i = x->size[0];
    x->size[0] = y->size[0];
    emxEnsureCapacity_real_T(x, i);
    nx = y->size[0] - 1;
    for (i = 0; i <= nx; i++) {
      x->data[i] = y->data[i];
    }

    emxInit_boolean_T(&r4, 1);
    b_mod(x, y);
    i = r4->size[0];
    r4->size[0] = y->size[0];
    emxEnsureCapacity_boolean_T(r4, i);
    nx = y->size[0];
    emxFree_real_T(&x);
    for (i = 0; i < nx; i++) {
      r4->data[i] = (y->data[i] == 0.0);
    }

    k = r4->size[0];
    for (i = 0; i < k; i++) {
      if (r4->data[i] && positiveInput->data[i]) {
        y->data[i] = 6.2831853071795862;
      }
    }

    emxFree_boolean_T(&r4);
    emxFree_boolean_T(&positiveInput);
    k = r->size[0];
    nx = 0;
    for (i = 0; i < k; i++) {
      if (r->data[i] || r1->data[i]) {
        xi->data[i] = y->data[nx] - 3.1415926535897931;
        nx++;
      }
    }

    nx = cosDelta->size[0];
    for (i = 0; i < nx; i++) {
      cosDelta->data[i] = az_tmp * cosDelta->data[i] + cosphi1_tmp *
        sinDelta->data[i] * cosAz;
    }

    nx = cosDelta->size[0];
    for (k = 0; k < nx; k++) {
      cosDelta->data[k] = asin(cosDelta->data[k]);
    }

    nx = cosDelta->size[0];
    for (i = 0; i < nx; i++) {
      cosDelta->data[i] *= 57.295779513082323;
    }

    nx = xi->size[0];
    for (i = 0; i < nx; i++) {
      xi->data[i] *= 57.295779513082323;
    }

    /*   Use exact endpoint. */
    cosDelta->data[cosDelta->size[0] - 1] = lat2;
    xi->data[xi->size[0] - 1] = lon2;
  } else {
    i = cosDelta->size[0];
    cosDelta->size[0] = 2;
    emxEnsureCapacity_real_T(cosDelta, i);
    cosDelta->data[0] = lat1;
    cosDelta->data[1] = lat2;
    i = xi->size[0];
    xi->size[0] = 2;
    emxEnsureCapacity_real_T(xi, i);
    xi->data[0] = lon1;
    xi->data[1] = lon2;
    w = 0.017453292519943295 * lat1;
    cosAz = 0.017453292519943295 * lat2;

    /* -------------------------------------------------------------------------- */
    /*  Great circle distance and azimuth between points on a sphere, using the */
    /*  Haversine Formula for distance.  All angles are in radians. */
    a = sin((cosAz - w) / 2.0);
    b_a = sin((0.017453292519943295 * lon2 - 0.017453292519943295 * lon1) / 2.0);
    i = dist->size[0];
    dist->size[0] = 1;
    emxEnsureCapacity_real_T(dist, i);
    dist->data[0] = 2.0 * asin(sqrt(a * a + cos(w) * cos(cosAz) * (b_a * b_a)));
  }

  i = lattrk->size[0];
  lattrk->size[0] = cosDelta->size[0];
  emxEnsureCapacity_real_T(lattrk, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    lattrk->data[i] = cosDelta->data[i];
  }

  i = lontrk->size[0];
  lontrk->size[0] = xi->size[0];
  emxEnsureCapacity_real_T(lontrk, i);
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    lontrk->data[i] = xi->data[i];
  }

  nx = dist->size[0];
  for (i = 0; i < nx; i++) {
    dist->data[i] *= actualradius;
  }

  /*  Elevation profile between the start and end points. */
  /* -------------------------------------------------------------------------- */
  /*  Use the griddentInterpolant object F which is defined in the intrinsic */
  /*  coordinate system referred to by raster reference object R. */
  /*  Interpolate bilinearly in intrinsic coordinates. */
  /* -------------------------------------------------------------------------- */
  /* latitudeToIntrinsicY Convert from latitude to intrinsic Y */
  /*  */
  /*    yIntrinsic = latitudeToIntrinsicY(R, LAT) returns the */
  /*    intrinsic Y value of the line corresponding to the small */
  /*    circle at latitude LAT, based on the relationship */
  /*    defined by the referencing object R. The input may */
  /*    include values that fall completely outside the latitude */
  /*    limits of the raster (or image). In this case yIntrinsic */
  /*    is either extrapolated outside the intrinsic Y limits -- */
  /*    for elements of LAT that fall within the interval */
  /*    [-90 90] degrees, or set to NaN -- for elements of LAT */
  /*    that do not correspond to valid latitudes. NaN-valued */
  /*    elements of LAT map to NaNs in yIntrinsic. */
  /*  Elements of LAT are less than -90 degrees or */
  /*  that exceed +90 degrees should map to NaN. */
  i = r->size[0];
  r->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r->data[i] = (cosDelta->data[i] < -90.0);
  }

  i = r1->size[0];
  r1->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r1, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r1->data[i] = (90.0 < cosDelta->data[i]);
  }

  k = r->size[0];
  for (i = 0; i < k; i++) {
    if (r->data[i] || r1->data[i]) {
      cosDelta->data[i] = rtNaN;
    }
  }

  /*  Shift and scale latitude */
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    cosDelta->data[i] = public_raster_struct->YIntrinsicLimits[0] +
      (cosDelta->data[i] - public_raster_struct->FirstCornerLatitude) *
      public_raster_struct->DeltaLatitudeDenominator /
      public_raster_struct->DeltaLatitudeNumerator;
  }

  /* -------------------------------------------------------------------------- */
  /* longitudeToIntrinsicX Convert from longitude to intrinsic X */
  /*  */
  /*    xIntrinsic = longitudeToIntrinsicX(R, LON) returns the */
  /*    intrinsic X value of the line corresponding to the */
  /*    meridian at longitude LON, based on the relationship */
  /*    defined by the referencing object R. The input may */
  /*    include values that fall completely outside the */
  /*    longitude limits of the raster (or image). In this case */
  /*    xIntrinsic is extrapolated outside the intrinsic X */
  /*    limits. NaN-valued elements of LON map to NaNs in */
  /*    xIntrinsic. */
  w = public_raster_struct->pLongitudeLimits[0];
  cosAz = public_raster_struct->pLongitudeLimits[1];

  /*  Adjust longitude wrapping to get within the limits, */
  /*  whenever possible. */
  if (public_raster_struct->pLongitudeLimits[1] -
      public_raster_struct->pLongitudeLimits[0] <=
      public_raster_struct->FullCycle) {
    if (public_raster_struct->DeltaLongitudeNumerator > 0.0) {
      /*  Wrap to interval R.FirstCornerLongitude + [0 360] */
      /* lon = w + R.WrapToCycleFcn(lon - w); etc */
      i = sinDelta->size[0];
      sinDelta->size[0] = xi->size[0];
      emxEnsureCapacity_real_T(sinDelta, i);
      nx = xi->size[0];
      for (i = 0; i < nx; i++) {
        sinDelta->data[i] = xi->data[i] - public_raster_struct->
          pLongitudeLimits[0];
      }

      c_mod(sinDelta, public_raster_struct->FullCycle, xi);
      nx = xi->size[0];
      for (i = 0; i < nx; i++) {
        xi->data[i] += public_raster_struct->pLongitudeLimits[0];
      }
    } else {
      /*  Wrap to interval R.FirstCornerLongitude + [-360 0] */
      i = sinDelta->size[0];
      sinDelta->size[0] = xi->size[0];
      emxEnsureCapacity_real_T(sinDelta, i);
      nx = xi->size[0];
      for (i = 0; i < nx; i++) {
        sinDelta->data[i] = public_raster_struct->pLongitudeLimits[1] - xi->
          data[i];
      }

      c_mod(sinDelta, public_raster_struct->FullCycle, xi);
      nx = xi->size[0];
      for (i = 0; i < nx; i++) {
        xi->data[i] = public_raster_struct->pLongitudeLimits[1] - xi->data[i];
      }
    }
  } else {
    /*  Any longitude can be wrapped to fall within the */
    /*  interval [w e], and in fact there's more than one */
    /*  solution for certain longitudes. Resolve the ambiguity */
    /*  by moving longitudes that are west of the western */
    /*  limit the minimal number of cycles to the east that */
    /*  puts them within the limits. Likewise, move longitudes */
    /*  that exceed the eastern limit the minimum number of */
    /*  cycles to the west. */
    k = xi->size[0] - 1;
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (xi->data[i] < w) {
        nx++;
      }
    }

    emxInit_int32_T(&r3, 1);
    i = r3->size[0];
    r3->size[0] = nx;
    emxEnsureCapacity_int32_T(r3, i);
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (xi->data[i] < w) {
        r3->data[nx] = i + 1;
        nx++;
      }
    }

    i = sinDelta->size[0];
    sinDelta->size[0] = r3->size[0];
    emxEnsureCapacity_real_T(sinDelta, i);
    nx = r3->size[0];
    for (i = 0; i < nx; i++) {
      sinDelta->data[i] = xi->data[r3->data[i] - 1] -
        public_raster_struct->pLongitudeLimits[0];
    }

    emxFree_int32_T(&r3);
    c_mod(sinDelta, public_raster_struct->FullCycle, y);
    k = xi->size[0];
    nx = 0;
    for (i = 0; i < k; i++) {
      if (xi->data[i] < w) {
        xi->data[i] = w + y->data[nx];
        nx++;
      }
    }

    w = public_raster_struct->pLongitudeLimits[1] -
      public_raster_struct->FullCycle;
    k = xi->size[0] - 1;
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (xi->data[i] > cosAz) {
        nx++;
      }
    }

    emxInit_int32_T(&r5, 1);
    i = r5->size[0];
    r5->size[0] = nx;
    emxEnsureCapacity_int32_T(r5, i);
    nx = 0;
    for (i = 0; i <= k; i++) {
      if (xi->data[i] > cosAz) {
        r5->data[nx] = i + 1;
        nx++;
      }
    }

    i = sinDelta->size[0];
    sinDelta->size[0] = r5->size[0];
    emxEnsureCapacity_real_T(sinDelta, i);
    nx = r5->size[0];
    for (i = 0; i < nx; i++) {
      sinDelta->data[i] = xi->data[r5->data[i] - 1] - w;
    }

    emxFree_int32_T(&r5);
    c_mod(sinDelta, public_raster_struct->FullCycle, y);
    k = xi->size[0];
    nx = 0;
    for (i = 0; i < k; i++) {
      if (xi->data[i] > cosAz) {
        xi->data[i] = w + y->data[nx];
        nx++;
      }
    }
  }

  emxFree_real_T(&y);
  emxFree_real_T(&sinDelta);

  /*  Shift and scale longitude */
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    xi->data[i] = public_raster_struct->XIntrinsicLimits[0] + (xi->data[i] -
      public_raster_struct->FirstCornerLongitude) *
      public_raster_struct->DeltaLongitudeDenominator /
      public_raster_struct->DeltaLongitudeNumerator;
  }

  /*  Snap in all points that fall within distance 0.5 of an edge, so that */
  /*  we get a non-NaN value for them from interp2. */
  i = r->size[0];
  r->size[0] = xi->size[0];
  emxEnsureCapacity_boolean_T(r, i);
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    r->data[i] = (0.5 <= xi->data[i]);
  }

  i = r1->size[0];
  r1->size[0] = xi->size[0];
  emxEnsureCapacity_boolean_T(r1, i);
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    r1->data[i] = (xi->data[i] < 1.0);
  }

  k = r->size[0];
  for (i = 0; i < k; i++) {
    if (r->data[i] && r1->data[i]) {
      xi->data[i] = 1.0;
    }
  }

  i = r->size[0];
  r->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r->data[i] = (0.5 <= cosDelta->data[i]);
  }

  i = r1->size[0];
  r1->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r1, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r1->data[i] = (cosDelta->data[i] < 1.0);
  }

  k = r->size[0];
  for (i = 0; i < k; i++) {
    if (r->data[i] && r1->data[i]) {
      cosDelta->data[i] = 1.0;
    }
  }

  i = r->size[0];
  r->size[0] = xi->size[0];
  emxEnsureCapacity_boolean_T(r, i);
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    r->data[i] = (public_raster_struct->RasterSize[1] < xi->data[i]);
  }

  i = r1->size[0];
  r1->size[0] = xi->size[0];
  emxEnsureCapacity_boolean_T(r1, i);
  nx = xi->size[0];
  for (i = 0; i < nx; i++) {
    r1->data[i] = (xi->data[i] <= public_raster_struct->RasterSize[1] + 0.5);
  }

  k = r->size[0];
  for (i = 0; i < k; i++) {
    if (r->data[i] && r1->data[i]) {
      xi->data[i] = public_raster_struct->RasterSize[1];
    }
  }

  i = r->size[0];
  r->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r->data[i] = (public_raster_struct->RasterSize[0] < cosDelta->data[i]);
  }

  i = r1->size[0];
  r1->size[0] = cosDelta->size[0];
  emxEnsureCapacity_boolean_T(r1, i);
  nx = cosDelta->size[0];
  for (i = 0; i < nx; i++) {
    r1->data[i] = (cosDelta->data[i] <= public_raster_struct->RasterSize[0] +
                   0.5);
  }

  k = r->size[0];
  for (i = 0; i < k; i++) {
    if (r->data[i] && r1->data[i]) {
      cosDelta->data[i] = public_raster_struct->RasterSize[0];
    }
  }

  emxFree_boolean_T(&r1);
  emxFree_boolean_T(&r);

  /*  testF = griddedInterpolant(Z); */
  /*  testF.ExtrapolationMethod = 'none'; */
  /*  v = testF(yi, xi); */
  /* linearly interpolates inside of elevation map Z */
  interp2_local(Z, xi, cosDelta, h);

  /*  Visibility of points along the profile between the start and end points. */
  /*  testF = griddedInterpolant(Z); */
  /*  testF.ExtrapolationMethod = 'none'; */
  /*  testv = testF(yi, xi); */
  /*   */
  /*  if testv ~= v */
  /*      fprintf("ow!"); */
  /*  end */
  /* -------------------------------------------------------------------------- */
  b = rtIsInf(apparentradius);
  emxFree_real_T(&cosDelta);
  emxFree_real_T(&xi);
  emxInit_real_T(&ang2, 2);
  emxInit_real_T(&b_dist, 2);
  if (!b) {
    i = b_dist->size[0] * b_dist->size[1];
    b_dist->size[0] = 1;
    b_dist->size[1] = dist->size[0];
    emxEnsureCapacity_real_T(b_dist, i);
    nx = dist->size[0];
    for (i = 0; i < nx; i++) {
      b_dist->data[i] = dist->data[i];
    }

    i = ang2->size[0] * ang2->size[1];
    ang2->size[0] = 1;
    ang2->size[1] = h->size[0];
    emxEnsureCapacity_real_T(ang2, i);
    nx = h->size[0];
    for (i = 0; i < nx; i++) {
      ang2->data[i] = h->data[i];
    }

    adjustterrain(b_dist, ang2, apparentradius, x1, z1);
  } else {
    i = x1->size[0] * x1->size[1];
    x1->size[0] = 1;
    x1->size[1] = dist->size[0];
    emxEnsureCapacity_real_T(x1, i);
    nx = dist->size[0];
    for (i = 0; i < nx; i++) {
      x1->data[i] = dist->data[i];
    }

    i = z1->size[0] * z1->size[1];
    z1->size[0] = 1;
    z1->size[1] = h->size[0];
    emxEnsureCapacity_real_T(z1, i);
    nx = h->size[0];
    for (i = 0; i < nx; i++) {
      z1->data[i] = h->data[i];
    }
  }

  /*  Convert AGL observer altitude to MSL  */
  /*   Observer is at first location */
  w = z1->data[0];

  /*  Shift terrain so observer is at altitude 0, and terrain altitudes are relative */
  /*  to the observer */
  cosAz = z1->data[0];
  nx = z1->size[0] * z1->size[1] - 1;
  i = z1->size[0] * z1->size[1];
  z1->size[0] = 1;
  emxEnsureCapacity_real_T(z1, i);
  for (i = 0; i <= nx; i++) {
    z1->data[i] -= cosAz;
  }

  /*  Observer is at first location */
  /*  Compute the angles of sight from the observer to each point on the profile. */
  /*  measured positive up from the center of the sphere */
  i = ang->size[0] * ang->size[1];
  ang->size[0] = 1;
  if (z1->size[1] <= x1->size[1]) {
    ang->size[1] = z1->size[1];
  } else {
    ang->size[1] = x1->size[1];
  }

  emxEnsureCapacity_real_T(ang, i);
  if (z1->size[1] <= x1->size[1]) {
    nx = z1->size[1];
  } else {
    nx = x1->size[1];
  }

  for (k = 0; k < nx; k++) {
    ang->data[k] = rt_atan2d_snf(z1->data[k], x1->data[k]);
  }

  i = ang->size[0] * ang->size[1];
  nx = ang->size[0] * ang->size[1];
  ang->size[0] = 1;
  emxEnsureCapacity_real_T(ang, nx);
  nx = i - 1;
  for (i = 0; i <= nx; i++) {
    ang->data[i] += 3.1415926535897931;
  }

  if ((x1->data[0] == 0.0) && (z1->data[0] == 0.0)) {
    ang->data[0] = 1.5707963267948966;

    /*  Look straight down at observer's location */
  }

  /*  Find the cumulative maximum:  maxtohere(k) equals max(ang(1:k)) */
  if ((ang->size[1] != 0) && (ang->size[1] != 1)) {
    i = ang->size[1];
    for (k = 0; k <= i - 2; k++) {
      ang->data[k + 1] = fmax(ang->data[k], ang->data[k + 1]);
    }
  }

  /*  Adjust the angles for the altitude of the target height above ground level  */
  /*  or sea level and redo calculation of angles. This makes the obscuring factor */
  /*  the terrain only, and not any target height. To model stuff above the terrain  */
  /*  like a forest canopy, pass in a z vector that has the added heights. */
  if (!b) {
    i = b_dist->size[0] * b_dist->size[1];
    b_dist->size[0] = 1;
    b_dist->size[1] = dist->size[0];
    emxEnsureCapacity_real_T(b_dist, i);
    nx = dist->size[0];
    for (i = 0; i < nx; i++) {
      b_dist->data[i] = dist->data[i];
    }

    i = ang2->size[0] * ang2->size[1];
    ang2->size[0] = 1;
    ang2->size[1] = h->size[0];
    emxEnsureCapacity_real_T(ang2, i);
    nx = h->size[0];
    for (i = 0; i < nx; i++) {
      ang2->data[i] = talt;
    }

    adjustterrain(b_dist, ang2, apparentradius, x2, z2);
    i = z2->size[0] * z2->size[1];
    nx = z2->size[0] * z2->size[1];
    z2->size[0] = 1;
    emxEnsureCapacity_real_T(z2, nx);
    nx = i - 1;
    for (i = 0; i <= nx; i++) {
      z2->data[i] -= w;
    }
  } else {
    i = z2->size[0] * z2->size[1];
    z2->size[0] = 1;
    z2->size[1] = h->size[0];
    emxEnsureCapacity_real_T(z2, i);
    nx = h->size[0];
    for (i = 0; i < nx; i++) {
      z2->data[i] = talt - w;
    }

    i = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = x1->size[1];
    emxEnsureCapacity_real_T(x2, i);
    nx = x1->size[0] * x1->size[1];
    for (i = 0; i < nx; i++) {
      x2->data[i] = x1->data[i];
    }
  }

  emxFree_real_T(&b_dist);

  /*  Compute line of sight angles again. */
  i = ang2->size[0] * ang2->size[1];
  ang2->size[0] = 1;
  if (z2->size[1] <= x2->size[1]) {
    ang2->size[1] = z2->size[1];
  } else {
    ang2->size[1] = x2->size[1];
  }

  emxEnsureCapacity_real_T(ang2, i);
  if (z2->size[1] <= x2->size[1]) {
    nx = z2->size[1];
  } else {
    nx = x2->size[1];
  }

  for (k = 0; k < nx; k++) {
    ang2->data[k] = rt_atan2d_snf(z2->data[k], x2->data[k]);
  }

  i = ang2->size[0] * ang2->size[1];
  nx = ang2->size[0] * ang2->size[1];
  ang2->size[0] = 1;
  emxEnsureCapacity_real_T(ang2, nx);
  nx = i - 1;
  for (i = 0; i <= nx; i++) {
    ang2->data[i] += 3.1415926535897931;
  }

  if ((x2->data[0] == 0.0) && (z2->data[0] == 0.0)) {
    ang2->data[0] = 1.5707963267948966;

    /*  Look straight down at observer's location */
  }

  emxInit_boolean_T(&vis, 2);

  /*  Visible are points that rise above the maximum angle of LOS of intervening  */
  /*  terrain. */
  i = vis->size[0] * vis->size[1];
  vis->size[0] = 1;
  vis->size[1] = ang2->size[1];
  emxEnsureCapacity_boolean_T(vis, i);
  nx = ang2->size[0] * ang2->size[1];
  for (i = 0; i < nx; i++) {
    vis->data[i] = (ang2->data[i] >= ang->data[i]);
  }

  emxFree_real_T(&ang2);
  emxFree_real_T(&ang);

  /*  Visibility of first point below terrain needs a special test, since */
  /*  it always passes the "angles of terrain up to here" test. */
  if ((z2->data[0] < z1->data[0]) && (z1->data[0] < 0.0)) {
    vis->data[0] = false;
  }

  i = vis0->size[0];
  vis0->size[0] = vis->size[1];
  emxEnsureCapacity_boolean_T(vis0, i);
  nx = vis->size[1];
  for (i = 0; i < nx; i++) {
    vis0->data[i] = vis->data[i];
  }

  emxFree_boolean_T(&vis);
}

/* End of code generation (mexable_calculateLOS.c) */
