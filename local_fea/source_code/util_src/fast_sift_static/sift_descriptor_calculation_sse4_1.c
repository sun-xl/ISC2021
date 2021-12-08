/*
 * sift_descriptor_calculation_sse4_1.c
 *
 *  Created on: Nov 10, 2014
 *      Author: peter
 */

#include "sift_descriptor_calculation.h"
#include "sift_math.h"
#include "sift_sse4_2_typedef.h"
#include <string.h>

void sift_calc_keypoint_descriptor_f_sse4_1(const SIFT_FILTER *f, float *descr,
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

  // speedup
  /*
  if (f->o_cur != k->o ||  f->o_cur != f->grad_o)
  {
          printf ("sift_calc_keypoint_descriptor_f_sse4_1: error\n");
          return;
  }
  */

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
  // Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0).
  float *const dpt = descr + SPATIAL_BINS2 * binyo + SPATIAL_BINS2 * binxo;

  const float st0 = sinf(angle0);
  const float ct0 = cosf(angle0);
  const float SBP = f->magnif * sigma + EPSILON_F;

  const float st0_SBP = st0 / SBP;
  const float ct0_SBP = ct0 / SBP;

  const float ORIENTATION_BINS_2PI = ORIENTATION_BINS / ((float)PI2);
  const float inverse_wsigma2_2_EXPN_SZ_EXPN_MAX =
      1.0 / (2.0 * f->windowSize * f->windowSize) * ((float)EXPN_SZ_EXPN_MAX);

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

  const __m128 onef_v128 = _mm_set1_ps(1.0f);
  const __m128i onei_v128 = _mm_set1_epi32(1);
  const __m128 halff_v128 = _mm_set1_ps(0.5f);
  const __m128i SPATIAL_BINS2_v128 = _mm_set1_epi32(SPATIAL_BINS2);
  const __m128i minus_SPATIAL_BINS2_v128 =
      _mm_set1_epi32(-SPATIAL_BINS2 - 1);  // be careful !!!
  const __m128 st0_SBP_v128 = _mm_set1_ps(st0_SBP);
  const __m128 ct0_SBP_v128 = _mm_set1_ps(ct0_SBP);
  const __m128 ORIENTATION_BINS_2PI_v128 = _mm_set1_ps(ORIENTATION_BINS_2PI);
  const __m128 inverse_wsigma2_2_EXPN_SZ_EXPN_MAX_v128 =
      _mm_set1_ps(inverse_wsigma2_2_EXPN_SZ_EXPN_MAX);
  const __m128 angle0_v128 = _mm_set1_ps(angle0);

  //	Process pixels in the intersection of the image rectangle
  //	(1,1)-(M-1,N-1) and the keypoint bounding box.
  for (dyi = SIFT_MAX(-W, 1 - yi);
       dyi <= SIFT_MIN(W, f->octave_height - yi - 2); dyi++) {
    const int lb = SIFT_MAX(-W, 1 - xi);
    const int hb = SIFT_MIN(W, f->octave_width - xi - 2);

    // y fractional displacement
    const float dy0 = yi_y + dyi;
    const __m128 dy128 = _mm_set1_ps(dy0);

    for (dxi = lb; dxi <= hb - 3; dxi += 4)  // sse 128
    {
      int p;
      const float *const mod_angle_ptr = pt + dxi * xo + dyi * yo;

      __m128 mod, angle, theta;
      __m128 data0, data1;
      // a1-m1-a0-m0
      // a3-m3-a2-m2
      // shuffle
      // shuffle pattern binary: 10001000
      // shuffle pattern binary: 11011101
      // m3-m2-m1-m0
      // a3-a2-a1-a0
      if (((uintptr_t)(mod_angle_ptr)) & 0xF)  // unaligned
      {
        data0 = _mm_loadu_ps(mod_angle_ptr);
        data1 = _mm_loadu_ps(mod_angle_ptr + 4);
      } else  // aligned
      {
        data0 = _mm_load_ps(mod_angle_ptr);
        data1 = _mm_load_ps(mod_angle_ptr + 4);
      }
      mod = _mm_shuffle_ps(data0, data1, 0x88);
      angle = _mm_shuffle_ps(data0, data1, 0xdd);

      theta = _mm_sub_ps(angle, angle0_v128);
      __m128 PI2_128 = _mm_set1_ps((float)PI2);
      __m128 ZERO_128 = _mm_setzero_ps();
      // clamping
      while (1) {
        __m128 cmp_res = _mm_cmpgt_ps(theta, PI2_128);
        int mask = _mm_movemask_ps(cmp_res);
        if (!mask) break;
        theta = _mm_sub_ps(theta, _mm_and_ps(cmp_res, PI2_128));
      }
      while (1) {
        __m128 cmp_res = _mm_cmplt_ps(theta, ZERO_128);
        int mask = _mm_movemask_ps(cmp_res);
        if (!mask) break;
        theta = _mm_add_ps(theta, _mm_and_ps(cmp_res, PI2_128));
      }

      // x fractional displacement
      const float dx0 = xi_x + dxi;
      const __m128 dx128 = _mm_set_ps(dx0 + 3, dx0 + 2, dx0 + 1, dx0);

      // get the displacement normalized (rotation angle0)
      // with regard to the keypoint orientation and extension
      __m128 nx, ny, nt;
      nx = _mm_add_ps(_mm_mul_ps(ct0_SBP_v128, dx128),
                      _mm_mul_ps(st0_SBP_v128, dy128));
      ny = _mm_sub_ps(_mm_mul_ps(ct0_SBP_v128, dy128),
                      _mm_mul_ps(st0_SBP_v128, dx128));
      nt = _mm_mul_ps(ORIENTATION_BINS_2PI_v128, theta);

      // Get the Gaussian weight of the sample. The Gaussian window
      // has a standard deviation equal to NBP/2. Note that dx and dy
      // are in the normalized frame, so that -NBP/2 <= dx <= NBP/2
      __m128 x = _mm_mul_ps(_mm_add_ps(_mm_mul_ps(nx, nx), _mm_mul_ps(ny, ny)),
                            inverse_wsigma2_2_EXPN_SZ_EXPN_MAX_v128);
      // floorf
      __m128 xfloor_f = _mm_floor_ps(x);
      __m128 r = _mm_sub_ps(x, xfloor_f);
      __m128 r1 = _mm_sub_ps(onef_v128, r);

      UV128I32 xfloor_i;
      xfloor_i.v = _mm_cvtps_epi32(xfloor_f);

      __m128 a = _mm_set_ps(
          f->fast_expn_tab[xfloor_i.i32[3]], f->fast_expn_tab[xfloor_i.i32[2]],
          f->fast_expn_tab[xfloor_i.i32[1]], f->fast_expn_tab[xfloor_i.i32[0]]);
      xfloor_i.v = _mm_add_epi32(xfloor_i.v, onei_v128);  // +1
      __m128 b = _mm_set_ps(
          f->fast_expn_tab[xfloor_i.i32[3]], f->fast_expn_tab[xfloor_i.i32[2]],
          f->fast_expn_tab[xfloor_i.i32[1]], f->fast_expn_tab[xfloor_i.i32[0]]);
      __m128 win = _mm_add_ps(_mm_mul_ps(r1, a), _mm_mul_ps(r, b));
      __m128 win_mod = _mm_mul_ps(win, mod);

      // The sample will be distributed in 8 adjacent bins.
      // We start from the ``lower-left'' bin.
      __m128 binxf, binyf, bintf;
      binxf = _mm_floor_ps(_mm_sub_ps(nx, halff_v128));
      binyf = _mm_floor_ps(_mm_sub_ps(ny, halff_v128));
      bintf = _mm_floor_ps(nt);

      __m128 rbinx, rbiny, rbint;
      rbinx = _mm_sub_ps(nx, _mm_add_ps(binxf, halff_v128));
      rbiny = _mm_sub_ps(ny, _mm_add_ps(binyf, halff_v128));
      rbint = _mm_sub_ps(nt, bintf);
      __m128 new_rbinx0, new_rbiny0, new_rbint0;
      new_rbinx0 = _mm_sub_ps(onef_v128, rbinx);
      new_rbiny0 = _mm_sub_ps(onef_v128, rbiny);
      new_rbint0 = _mm_sub_ps(onef_v128, rbint);

      UV128I32 binx, biny, bint;
      binx.v = _mm_cvtps_epi32(binxf);
      biny.v = _mm_cvtps_epi32(binyf);
      bint.v = _mm_cvtps_epi32(bintf);
      UV128I32 new_binx1, new_biny1, new_bint1;
      new_binx1.v = _mm_add_epi32(binx.v, onei_v128);
      new_biny1.v = _mm_add_epi32(biny.v, onei_v128);
      new_bint1.v = _mm_add_epi32(bint.v, onei_v128);

      __m128i binx_flag, biny_flag, new_binx1_flag, new_biny1_flag;
      binx_flag =
          _mm_and_si128(_mm_cmpgt_epi32(binx.v, minus_SPATIAL_BINS2_v128),
                        _mm_cmplt_epi32(binx.v, SPATIAL_BINS2_v128));
      biny_flag =
          _mm_and_si128(_mm_cmpgt_epi32(biny.v, minus_SPATIAL_BINS2_v128),
                        _mm_cmplt_epi32(biny.v, SPATIAL_BINS2_v128));
      new_binx1_flag =
          _mm_and_si128(_mm_cmpgt_epi32(new_binx1.v, minus_SPATIAL_BINS2_v128),
                        _mm_cmplt_epi32(new_binx1.v, SPATIAL_BINS2_v128));
      new_biny1_flag =
          _mm_and_si128(_mm_cmpgt_epi32(new_biny1.v, minus_SPATIAL_BINS2_v128),
                        _mm_cmplt_epi32(new_biny1.v, SPATIAL_BINS2_v128));

      UV128I32 binx_biny_flag, new_binx1_biny_flag, binx_new_bin1_flag,
          new_binx1_new_biny1_flag;
      binx_biny_flag.v = _mm_and_si128(binx_flag, biny_flag);
      new_binx1_biny_flag.v = _mm_and_si128(new_binx1_flag, biny_flag);
      binx_new_bin1_flag.v = _mm_and_si128(binx_flag, new_biny1_flag);
      new_binx1_new_biny1_flag.v =
          _mm_and_si128(new_binx1_flag, new_biny1_flag);

      // Distribute the current sample into the 8 adjacent bins
      if (_mm_movemask_epi8(binx_biny_flag.v)) {
        __m128 weight = _mm_mul_ps(_mm_mul_ps(win_mod, new_rbinx0), new_rbiny0);
        UV128F32 weight0, weight1;
        weight0.v = _mm_mul_ps(weight, new_rbint0);
        weight1.v = _mm_mul_ps(weight, rbint);

// t=0, x=0, y=0
// t=1, x=0, y=0
#define sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1(p)                \
  if (binx_biny_flag.i32[p] != 0) {                                       \
    atd(binx.i32[p], biny.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                   \
    atd(binx.i32[p], biny.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                   \
  }
        sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1(0)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1(1)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1(2)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1(3)
#undef sift_calc_keypoint_descriptor_f_sse4_1_hist_0_1
      }

      if (_mm_movemask_epi8(new_binx1_biny_flag.v)) {
        __m128 weight = _mm_mul_ps(_mm_mul_ps(win_mod, rbinx), new_rbiny0);
        UV128F32 weight0, weight1;
        weight0.v = _mm_mul_ps(weight, new_rbint0);
        weight1.v = _mm_mul_ps(weight, rbint);

// t=0, x=1, y=0
// t=1, x=1, y=0
#define sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3(p)                     \
  if (new_binx1_biny_flag.i32[p] != 0) {                                       \
    atd(new_binx1.i32[p], biny.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                        \
    atd(new_binx1.i32[p], biny.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                        \
  }
        sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3(0)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3(1)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3(2)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3(3)
#undef sift_calc_keypoint_descriptor_f_sse4_1_hist_2_3
      }

      if (_mm_movemask_epi8(binx_new_bin1_flag.v)) {
        __m128 weight = _mm_mul_ps(_mm_mul_ps(win_mod, new_rbinx0), rbiny);
        UV128F32 weight0, weight1;
        weight0.v = _mm_mul_ps(weight, new_rbint0);
        weight1.v = _mm_mul_ps(weight, rbint);

// t=0, x=0, y=1
// t=1, x=0, y=1
#define sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5(p)                     \
  if (binx_new_bin1_flag.i32[p] != 0) {                                        \
    atd(binx.i32[p], new_biny1.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                        \
    atd(binx.i32[p], new_biny1.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                        \
  }
        sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5(0)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5(1)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5(2)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5(3)
#undef sift_calc_keypoint_descriptor_f_sse4_1_hist_4_5
      }

      if (_mm_movemask_epi8(new_binx1_new_biny1_flag.v)) {
        __m128 weight = _mm_mul_ps(_mm_mul_ps(win_mod, rbinx), rbiny);
        UV128F32 weight0, weight1;
        weight0.v = _mm_mul_ps(weight, new_rbint0);
        weight1.v = _mm_mul_ps(weight, rbint);
// t=0, x=1, y=1
// t=1, x=1, y=1
#define sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7(p)                     \
  if (new_binx1_new_biny1_flag.i32[p] != 0) {                                  \
    atd(new_binx1.i32[p], new_biny1.i32[p], bint.i32[p] % ORIENTATION_BINS) += \
        weight0.f32[p];                                                        \
    atd(new_binx1.i32[p], new_biny1.i32[p],                                    \
        new_bint1.i32[p] % ORIENTATION_BINS) += weight1.f32[p];                \
  }
        sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7(0)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7(1)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7(2)
        sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7(3)
#undef sift_calc_keypoint_descriptor_f_sse4_1_hist_6_7
      }
    }

    for (; dxi <= hb; dxi++)  // c
    {
      // retrieve
      float mod = *(pt + dxi * xo + dyi * yo);        // magnitude
      float angle = *(pt + dxi * xo + dyi * yo + 1);  // orientation
      float theta = angle - angle0;

      while (theta > (float)PI2) theta -= (float)PI2;
      while (theta < 0.0F) theta += (float)PI2;

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
      float x = (nx * nx + ny * ny) * inverse_wsigma2_2_EXPN_SZ_EXPN_MAX;
      float a, b, r;
      int i;

      // floor
      int xi = (int)x;
      if (x >= 0 || (float)xi == x)
        i = xi;
      else
        i = xi - 1;

      r = x - i;
      a = f->fast_expn_tab[i];
      b = f->fast_expn_tab[i + 1];
      float win = a + r * (b - a);  // (1-r)*a + r*b
      const float win_mod = win * mod;

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
        weight = win_mod * new_rbinx0 * new_rbiny0;
        atd(binx, biny, bint % ORIENTATION_BINS) += weight * new_rbint0;
        atd(binx, biny, new_bint1 % ORIENTATION_BINS) += weight * rbint;
      }

      // t=0, x=1, y=0
      // t=1, x=1, y=0
      if (new_binx1 >= -SPATIAL_BINS2 && new_binx1 < SPATIAL_BINS2 &&
          biny >= -SPATIAL_BINS2 && biny < SPATIAL_BINS2) {
        weight = win_mod * rbinx * new_rbiny0;
        atd(new_binx1, biny, bint % ORIENTATION_BINS) += weight * new_rbint0;
        atd(new_binx1, biny, new_bint1 % ORIENTATION_BINS) += weight * rbint;
      }

      // t=0, x=0, y=1
      // t=1, x=0, y=1
      if (binx >= -SPATIAL_BINS2 && binx < SPATIAL_BINS2 &&
          new_biny1 >= -SPATIAL_BINS2 && new_biny1 < SPATIAL_BINS2) {
        weight = win_mod * new_rbinx0 * rbiny;
        atd(binx, new_biny1, bint % ORIENTATION_BINS) += weight * new_rbint0;
        atd(binx, new_biny1, new_bint1 % ORIENTATION_BINS) += weight * rbint;
      }

      // t=0, x=1, y=1
      // t=1, x=1, y=1
      if (new_binx1 >= -SPATIAL_BINS2 && new_binx1 < SPATIAL_BINS2 &&
          new_biny1 >= -SPATIAL_BINS2 && new_biny1 < SPATIAL_BINS2) {
        weight = win_mod * rbinx * rbiny;
        atd(new_binx1, new_biny1, bint % ORIENTATION_BINS) +=
            weight * new_rbint0;
        atd(new_binx1, new_biny1, new_bint1 % ORIENTATION_BINS) +=
            weight * rbint;
      }
    }  // for (; dxi <= hb; dxi++)
  }    // for (dyi = MAX(-W, 1-yi); dyi<=MIN(W, f->octave_height-yi-2); dyi++)

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
