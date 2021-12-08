/*
 * sift_descriptor_calculation.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_descriptor_calculation.h"
#include "sift_math.h"
#include <string.h>

void sift_calc_keypoint_descriptor_f(const SIFT_FILTER *f, float *descr,
                                     const SIFT_KEYPOINT *k,
                                     const float angle0) {

  //    The SIFT descriptor is a three dimensional histogram of the
  //    position and orientation of the gradient.  There are SPATIAL_BINS bins
  // for
  //    each spatial dimension and ORIENTATION_BINS bins for the orientation
  // dimension,
  //    for a total of SPATIAL_BINS x SPATIAL_BINS x ORIENTATION_BINS bins.
  //
  //    The support of each spatial bin has an extension of SBP = 3sigma
  //    pixels, where sigma is the scale of the keypoint.  Thus all the
  //    bins together have a support SBP x NBP pixels wide. Since
  //    weighting and interpolation of pixel is used, the support extends
  //    by another half bin. Therefore, the support is a square window of
  //    SBP x (NBP + 1) pixels. Finally, since the patch can be
  //    arbitrarily rotated, we need to consider a window 2W += sqrt(2) x
  //    SBP x (NBP + 1) pixels wide.

  const int xo = 2;
  const int yo = 2 * f->octave_width;
  const int so = 2 * f->octave_size;

  const float x = k->x / f->xper;
  const float y = k->y / f->xper;
  const float sigma = k->sigma / f->xper;

  const int xi = (int)(x + 0.5);
  const int yi = (int)(y + 0.5);
  const int si = k->is;

  const float xi_x = xi - x;
  const float yi_y = yi - y;

  // center the scale space on the current keypoint
  const float *const pt =
      f->grad + xo * xi + yo * yi + so * (si - f->s_min - 1);

  const int binto = 1;
  const int binyo = SPATIAL_BINS * ORIENTATION_BINS;
  const int binxo = ORIENTATION_BINS;
  // center the descriptor on the current keypoint
  // Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0). // Peter
  // (NBP/2,NBP/2,0)
  float *const dpt =
      descr + (SPATIAL_BINS / 2) * binyo + (SPATIAL_BINS / 2) * binxo;

  const float st0 = sinf(angle0);
  const float ct0 = cosf(angle0);
  const float SBP = f->magnif * sigma + EPSILON_F;

  const float st0_SBP = st0 / SBP;
  const float ct0_SBP = ct0 / SBP;
  const float ORIENTATION_BINS_2PI = ORIENTATION_BINS / ((float)PI2);
  const float inverse_wsigma2_2 = 1.0 / (2.0 * f->windowSize * f->windowSize);

  const int W = (int)floorf(sqrtf(2.0) * SBP * (SPATIAL_BINS + 1) / 2.0 +
                            0.5);  // floor(3.0 * (4 + 1) / 2.0 + 0.5) = 8.0
  int dxi, dyi;

  // skip the keypoint if it is out of bounds
  // Peter: why the boundary is different???
  //	if (xi < 0 || xi > f->octave_width - 1
  //			|| yi < 0 || yi > f->octave_height - 1
  //			|| si < f->s_min + 1 || si > f->s_max - 2)

  if (xi < 0 || xi >= f->octave_width || yi < 0 || yi >= f->octave_height - 1 ||
      si < f->s_min + 1 || si > f->s_max - 2)
    return;

  // clear descriptor
  memset(descr, 0, sizeof(float) * DESCRIPTOR_SIZE);

#undef atd
#define atd(dbinx, dbiny, dbint) \
  *(dpt + (dbint) * binto + (dbiny) * binyo + (dbinx) * binxo)

  //	Process pixels in the intersection of the image rectangle
  //	(1,1)-(M-1,N-1) and the keypoint bounding box.
  for (dyi = SIFT_MAX(-W, 1 - yi);
       dyi <= SIFT_MIN(W, f->octave_height - yi - 2); dyi++) {
    for (dxi = SIFT_MAX(-W, 1 - xi);
         dxi <= SIFT_MIN(W, f->octave_width - xi - 2); dxi++) {
      // retrieve
      float mod = *(pt + dxi * xo + dyi * yo);        // magnitude
      float angle = *(pt + dxi * xo + dyi * yo + 1);  // orientation
      float theta = angle - angle0;

      while (theta > (float)(PI2)) theta -= (float)(PI2);
      while (theta < 0.0F) theta += (float)(PI2);

      // fractional displacement
      float dx = xi_x + dxi;
      float dy = yi_y + dyi;

      // get the displacement normalized (rotation angle0)
      // with regard to the keypoint orientation and extension
      float nx = ct0_SBP * dx + st0_SBP * dy;
      float ny = -st0_SBP * dx + ct0_SBP * dy;
      float nt = ORIENTATION_BINS_2PI * theta;

      // Get the Gaussian weight of the sample. The Gaussian window
      // has a standard deviation equal to NBP/2. Note that dx and dy
      // are in the normalized frame, so that -NBP/2 <= dx <= NBP/2
      float x = (nx * nx + ny * ny) * inverse_wsigma2_2;
      float a, b, r;
      int i;

      x *= (float)EXPN_SZ_EXPN_MAX;

      int xi = (int)x;
      if (x >= 0 || (float)xi == x)
        i = xi;
      else
        i = xi - 1;

      r = x - i;
      a = f->fast_expn_tab[i];
      b = f->fast_expn_tab[i + 1];
      float win = a + r * (b - a);

      // The sample will be distributed in 8 adjacent bins.
      // We start from the ``lower-left'' bin.
      int binx, biny, bint;
      float fval;
      int ival;

      fval = nx - 0.5;
      ival = (int)fval;
      if (fval >= 0 || (float)ival == fval)
        binx = ival;
      else
        binx = ival - 1;

      fval = ny - 0.5;
      ival = (int)fval;
      if (fval >= 0 || (float)ival == fval)
        biny = ival;
      else
        biny = ival - 1;

      fval = nt;
      ival = (int)fval;
      if (fval >= 0 || (float)ival == fval)
        bint = ival;
      else
        bint = ival - 1;

      float rbinx = nx - (binx + 0.5);
      float rbiny = ny - (biny + 0.5);
      float rbint = nt - bint;

      const float win_mod = win * mod;
      const int new_binx1 = binx + 1;
      const int new_biny1 = biny + 1;
      const int new_bint1 = bint + 1;

      const float new_rbinx0 = (1 - rbinx);
      const float new_rbiny0 = (1 - rbiny);
      const float new_rbint0 = (1 - rbint);

      float weight;
      // Distribute the current sample into the 8 adjacent bins
      // t=0, x=0, y=0
      // t=1, x=0, y=0
      if (binx >= -SPATIAL_BINS2 && binx < SPATIAL_BINS2 &&
          biny >= -SPATIAL_BINS2 && biny < SPATIAL_BINS2) {
        weight = win_mod * new_rbinx0 * new_rbiny0 * new_rbint0;
        atd(binx, biny, bint % ORIENTATION_BINS) += weight;

        weight = win_mod * new_rbinx0 * new_rbiny0 * rbint;
        atd(binx, biny, new_bint1 % ORIENTATION_BINS) += weight;
      }

      // t=0, x=1, y=0
      // t=1, x=1, y=0
      if (new_binx1 >= -SPATIAL_BINS2 && new_binx1 < SPATIAL_BINS2 &&
          biny >= -SPATIAL_BINS2 && biny < SPATIAL_BINS2) {
        weight = win_mod * rbinx * new_rbiny0 * new_rbint0;
        atd(new_binx1, biny, bint % ORIENTATION_BINS) += weight;

        weight = win_mod * rbinx * new_rbiny0 * rbint;
        atd(new_binx1, biny, new_bint1 % ORIENTATION_BINS) += weight;
      }

      // t=0, x=0, y=1
      // t=1, x=0, y=1
      if (binx >= -SPATIAL_BINS2 && binx < SPATIAL_BINS2 &&
          new_biny1 >= -SPATIAL_BINS2 && new_biny1 < SPATIAL_BINS2) {
        weight = win_mod * new_rbinx0 * rbiny * new_rbint0;
        atd(binx, new_biny1, bint % ORIENTATION_BINS) += weight;

        weight = win_mod * new_rbinx0 * rbiny * rbint;
        atd(binx, new_biny1, new_bint1 % ORIENTATION_BINS) += weight;
      }

      // t=0, x=1, y=1
      // t=1, x=1, y=1
      if (new_binx1 >= -SPATIAL_BINS2 && new_binx1 < SPATIAL_BINS2 &&
          new_biny1 >= -SPATIAL_BINS2 && new_biny1 < SPATIAL_BINS2) {
        weight = win_mod * rbinx * rbiny * new_rbint0;
        atd(new_binx1, new_biny1, bint % ORIENTATION_BINS) += weight;

        weight = win_mod * rbinx * rbiny * rbint;
        atd(new_binx1, new_biny1, new_bint1 % ORIENTATION_BINS) += weight;
      }
    }  // for (dxi=MAX(-W, 1-xi); dxi<=MIN(W, f->octave_width-xi-2); dxi++)

  }  // for (dyi = MAX(-W, 1-yi); dyi<=MIN(W, f->octave_height-yi-2); dyi++)

  // Standard SIFT descriptors are normalized, truncated and normalized again
  {
    int bin;
    // Normalize the histogram to L2 unit length.
    float norm = 0.0;
    for (bin = 0; bin < DESCRIPTOR_SIZE; bin++)
      norm += descr[bin] * descr[bin];  // L2 norm
    if (norm < 1e-8)
      norm = EPSILON_F;
    else {
      float x = norm;
      // 32-bit version
      union {
        float x;
        int i;
      } u;
      float xhalf = (float)0.5 * x;
      // convert floating point value in RAW integer
      u.x = x;
      // gives initial guess y0
      u.i = 0x5f3759df - (u.i >> 1);
      /* two Newton steps */
      u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
      u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
      norm *= u.x;
      norm += EPSILON_F;
    }
    for (bin = 0; bin < DESCRIPTOR_SIZE; bin++) descr[bin] /= norm;

    // Set the descriptor to zero if it is lower than our norm_threshold
    if (f->norm_thresh && norm < f->norm_thresh) {
      for (bin = 0; bin < DESCRIPTOR_SIZE; ++bin) descr[bin] = 0;
    } else {
      // Truncate at 0.2
      for (bin = 0; bin < DESCRIPTOR_SIZE; ++bin) {
        if (descr[bin] > 0.2) descr[bin] = 0.2;
      }

      // Normalize again
      float norm = 0.0;
      for (bin = 0; bin < DESCRIPTOR_SIZE; bin++)
        norm += descr[bin] * descr[bin];  // L2 norm
      if (norm < 1e-8)
        norm = EPSILON_F;
      else {
        float x = norm;
        // 32-bit version
        union {
          float x;
          int i;
        } u;
        float xhalf = (float)0.5 * x;
        // convert floating point value in RAW integer
        u.x = x;
        // gives initial guess y0
        u.i = 0x5f3759df - (u.i >> 1);
        /* two Newton steps */
        u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
        u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
        norm *= u.x;
        norm += EPSILON_F;
      }
      for (bin = 0; bin < DESCRIPTOR_SIZE; bin++) descr[bin] /= norm;
    }
  }

#undef atd
}
