/*
 * sift_descriptor_calculation_avx.c
 *
 *  Created on: Nov 10, 2014
 *      Author: peter
 */

#include "sift_descriptor_calculation.h"
#include "sift_math.h"
#include "sift_avx_typedef.h"
#include <string.h>

void sift_calc_keypoint_descriptor_f_avx(const SIFT_FILTER *f, float *descr,
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
  // Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0). // Peter
  // (NBP/2,NBP/2,0)
  float *const dpt = descr + SPATIAL_BINS2 * binyo + SPATIAL_BINS2 * binxo;

  const float st0 = sinf(angle0);
  const float ct0 = cosf(angle0);
  const float SBP = f->magnif * sigma + EPSILON_F;

  const float st0_SBP = st0 / SBP;
  const float ct0_SBP = ct0 / SBP;

  const float ORIENTATION_BINS_2PI = (float)ORIENTATION_BINS / (float)PI2;
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

  //	Process pixels in the intersection of the image rectangle
  //	(1,1)-(M-1,N-1) and the keypoint bounding box.
  // some constant
  const __m256 onef_v256 = _mm256_set1_ps(1.0f);
  const __m256 halff_v256 = _mm256_set1_ps(0.5f);
  // const __m256i onei_v256 = _mm256_set1_epi32 (1); // need, if avx2 is
  // supported
  const __m256 SPATIAL_BINS2_v256 = _mm256_set1_ps(SPATIAL_BINS2);
  const __m256 minus_SPATIAL_BINS2_v256 = _mm256_set1_ps(-SPATIAL_BINS2);
  const __m256 angle0_v256 = _mm256_set1_ps(angle0);
  const __m256 st0_SBP_v256 = _mm256_set1_ps(st0_SBP);
  const __m256 ct0_SBP_v256 = _mm256_set1_ps(ct0_SBP);
  const __m256 ORIENTATION_BINS_2PI_v256 = _mm256_set1_ps(ORIENTATION_BINS_2PI);
  const __m256 inverse_wsigma2_2_EXPN_SZ_EXPN_MAX_v256 =
      _mm256_set1_ps(inverse_wsigma2_2_EXPN_SZ_EXPN_MAX);

  for (dyi = SIFT_MAX(-W, 1 - yi);
       dyi <= SIFT_MIN(W, f->octave_height - yi - 2); dyi++) {
    const int lb = SIFT_MAX(-W, 1 - xi);
    const int hb = SIFT_MIN(W, f->octave_width - xi - 2);

    // y fractional displacement
    const float dy0 = yi_y + dyi;
    const __m256 dy = _mm256_set1_ps(dy0);

    for (dxi = lb; dxi <= hb - 7; dxi += 8)  // avx 256
    {
      int p;
      const float *const mod_angle_ptr = pt + dxi * xo + dyi * yo;
      __m256 mod, angle;

      // a3-m3-a2-m2-a1-m1-a0-m0
      // a7-m7-a6-m6-a5-m5-a4-m4
      // shuffle
      // binary shuffle pattern: 10001000 0x88
      // binary shuffle pattern: 11011101 0xdd
      // mod: m7-m6-m3-m2-m5-m4-m1-m0
      // angle: a7-a6-a3-a2-a5-a4-a1-a0
      if (((uintptr_t)(mod_angle_ptr)) & 0x1F)  // unaligned load
      {
        __m256 data0, data1;
        data0 = _mm256_loadu_ps(mod_angle_ptr);
        data1 = _mm256_loadu_ps(mod_angle_ptr + 8);
        mod = _mm256_shuffle_ps(data0, data1, 0x88);
        angle = _mm256_shuffle_ps(data0, data1, 0xdd);
      } else  // aligned load
      {
        __m256 data0, data1;
        data0 = _mm256_load_ps(mod_angle_ptr);
        data1 = _mm256_load_ps(mod_angle_ptr + 8);
        mod = _mm256_shuffle_ps(data0, data1, 0x88);
        angle = _mm256_shuffle_ps(data0, data1, 0xdd);
      }

      __m256 theta = _mm256_sub_ps(angle, angle0_v256);
      // clamping
      __m256 PI2_256 = _mm256_set1_ps((float)PI2);
      __m256 ZERO_256 = _mm256_setzero_ps();
      while (1) {
        __m256 cmp_res = _mm256_cmp_ps(theta, PI2_256, 14);  // gt
        int mask = _mm256_movemask_ps(cmp_res);
        if (!mask) break;
        theta = _mm256_sub_ps(theta, _mm256_and_ps(cmp_res, PI2_256));
      }

      while (1) {
        __m256 cmp_res = _mm256_cmp_ps(theta, ZERO_256, 1);  // lt
        int mask = _mm256_movemask_ps(cmp_res);
        if (!mask) break;
        theta = _mm256_add_ps(theta, _mm256_and_ps(cmp_res, PI2_256));
      }

      // x fractional displacement
      const float dx0 = xi_x + dxi;
      // be careful about the order: a7-a6-a3-a2-a5-a4-a1-a0
      const __m256 dx = _mm256_set_ps(dx0 + 7, dx0 + 6, dx0 + 3, dx0 + 2,
                                      dx0 + 5, dx0 + 4, dx0 + 1, dx0);

      // get the displacement normalized (rotation angle0)
      // with regard to the keypoint orientation and extension
      __m256 nx, ny, nt;
      nx = _mm256_add_ps(_mm256_mul_ps(ct0_SBP_v256, dx),
                         _mm256_mul_ps(st0_SBP_v256, dy));
      ny = _mm256_sub_ps(_mm256_mul_ps(ct0_SBP_v256, dy),
                         _mm256_mul_ps(st0_SBP_v256, dx));
      nt = _mm256_mul_ps(ORIENTATION_BINS_2PI_v256, theta);

      // Get the Gaussian weight of the sample. The Gaussian window
      // has a standard deviation equal to NBP/2. Note that dx and dy
      // are in the normalized frame, so that -NBP/2 <= dx <= NBP/2
      __m256 x = _mm256_mul_ps(
          _mm256_add_ps(_mm256_mul_ps(nx, nx), _mm256_mul_ps(ny, ny)),
          inverse_wsigma2_2_EXPN_SZ_EXPN_MAX_v256);
      __m256 xfloor_f = _mm256_floor_ps(x);
      __m256 r = _mm256_sub_ps(x, xfloor_f);
      __m256 r1 = _mm256_sub_ps(onef_v256, r);

      UV256I32 xfloor_i, xfloor_1_i;
      xfloor_i.v = _mm256_cvtps_epi32(xfloor_f);
      // use _mm256_add_epi32, if avx2 is supproted
      xfloor_1_i.v =
          _mm256_cvtps_epi32(_mm256_add_ps(onef_v256, xfloor_f));  // 6 cycles

      __m256 a = _mm256_set_ps(
          f->fast_expn_tab[xfloor_i.i32[7]], f->fast_expn_tab[xfloor_i.i32[6]],
          f->fast_expn_tab[xfloor_i.i32[5]], f->fast_expn_tab[xfloor_i.i32[4]],
          f->fast_expn_tab[xfloor_i.i32[3]], f->fast_expn_tab[xfloor_i.i32[2]],
          f->fast_expn_tab[xfloor_i.i32[1]], f->fast_expn_tab[xfloor_i.i32[0]]);
      __m256 b = _mm256_set_ps(f->fast_expn_tab[xfloor_1_i.i32[7]],
                               f->fast_expn_tab[xfloor_1_i.i32[6]],
                               f->fast_expn_tab[xfloor_1_i.i32[5]],
                               f->fast_expn_tab[xfloor_1_i.i32[4]],
                               f->fast_expn_tab[xfloor_1_i.i32[3]],
                               f->fast_expn_tab[xfloor_1_i.i32[2]],
                               f->fast_expn_tab[xfloor_1_i.i32[1]],
                               f->fast_expn_tab[xfloor_1_i.i32[0]]);
      __m256 win = _mm256_add_ps(_mm256_mul_ps(r1, a), _mm256_mul_ps(r, b));

      __m256 win_mod = _mm256_mul_ps(win, mod);

      __m256 binx_f, biny_f, bint_f;
      binx_f = _mm256_floor_ps(_mm256_sub_ps(nx, halff_v256));
      biny_f = _mm256_floor_ps(_mm256_sub_ps(ny, halff_v256));
      bint_f = _mm256_floor_ps(nt);
      __m256 new_binx1_f, new_biny1_f;
      new_binx1_f = _mm256_add_ps(binx_f, onef_v256);
      new_biny1_f = _mm256_add_ps(biny_f, onef_v256);

      // use _mm256_and_si256, if avx2 is supported
      __m256 binx_flag_f, biny_flag_f, new_binx1_flag_f, new_biny1_flag_f;
      binx_flag_f =
          _mm256_and_ps(_mm256_cmp_ps(binx_f, minus_SPATIAL_BINS2_v256, 13),
                        _mm256_cmp_ps(binx_f, SPATIAL_BINS2_v256, 1));
      biny_flag_f =
          _mm256_and_ps(_mm256_cmp_ps(biny_f, minus_SPATIAL_BINS2_v256, 13),
                        _mm256_cmp_ps(biny_f, SPATIAL_BINS2_v256, 1));
      new_binx1_flag_f = _mm256_and_ps(
          _mm256_cmp_ps(new_binx1_f, minus_SPATIAL_BINS2_v256, 13),
          _mm256_cmp_ps(new_binx1_f, SPATIAL_BINS2_v256, 1));
      new_biny1_flag_f = _mm256_and_ps(
          _mm256_cmp_ps(new_biny1_f, minus_SPATIAL_BINS2_v256, 13),
          _mm256_cmp_ps(new_biny1_f, SPATIAL_BINS2_v256, 1));

      UV256F32 binx_biny_flag_f, new_binx1_biny_flag_f, binx_new_biny1_flag_f,
          new_binx1_new_biny1_flag_f;
      binx_biny_flag_f.v = _mm256_and_ps(binx_flag_f, biny_flag_f);
      new_binx1_biny_flag_f.v = _mm256_and_ps(new_binx1_flag_f, biny_flag_f);
      binx_new_biny1_flag_f.v = _mm256_and_ps(binx_flag_f, new_biny1_flag_f);
      new_binx1_new_biny1_flag_f.v =
          _mm256_and_ps(new_binx1_flag_f, new_biny1_flag_f);

      UV256I32 binx, biny, bint;
      binx.v = _mm256_cvtps_epi32(binx_f);
      biny.v = _mm256_cvtps_epi32(biny_f);
      bint.v = _mm256_cvtps_epi32(bint_f);
      UV256I32 new_binx1, new_biny1, new_bint1;
      // use _mm256_add_epi32, if avx2 is supported
      new_binx1.v = _mm256_cvtps_epi32(_mm256_add_ps(binx_f, onef_v256));
      new_biny1.v = _mm256_cvtps_epi32(_mm256_add_ps(biny_f, onef_v256));
      new_bint1.v = _mm256_cvtps_epi32(_mm256_add_ps(bint_f, onef_v256));

      __m256 rbinx, rbiny, rbint;
      rbinx = _mm256_sub_ps(nx, _mm256_add_ps(binx_f, halff_v256));
      rbiny = _mm256_sub_ps(ny, _mm256_add_ps(biny_f, halff_v256));
      rbint = _mm256_sub_ps(nt, bint_f);

      __m256 new_rbinx0, new_rbiny0, new_rbint0;
      new_rbinx0 = _mm256_sub_ps(onef_v256, rbinx);
      new_rbiny0 = _mm256_sub_ps(onef_v256, rbiny);
      new_rbint0 = _mm256_sub_ps(onef_v256, rbint);

      // Distribute the current sample into the 8 adjacent bins
      if (_mm256_movemask_ps(binx_biny_flag_f.v)) {
        __m256 weight =
            _mm256_mul_ps(_mm256_mul_ps(win_mod, new_rbinx0), new_rbiny0);
        UV256F32 weight0, weight1;
        weight0.v = _mm256_mul_ps(weight, new_rbint0);
        weight1.v = _mm256_mul_ps(weight, rbint);

// t=0, x=0, y=0
// t=1, x=0, y=0
#define sift_calc_keypoint_descriptor_f_avx_hist_0_1(p)                   \
  if (binx_biny_flag_f.f32[p]) {                                          \
    atd(binx.i32[p], biny.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                   \
    atd(binx.i32[p], biny.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                   \
  }
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(0)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(1)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(2)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(3)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(4)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(5)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(6)
        sift_calc_keypoint_descriptor_f_avx_hist_0_1(7)
#undef sift_calc_keypoint_descriptor_f_avx_hist_0_1
      }

      if (_mm256_movemask_ps(new_binx1_biny_flag_f.v)) {
        __m256 weight =
            _mm256_mul_ps(_mm256_mul_ps(win_mod, rbinx), new_rbiny0);
        UV256F32 weight0, weight1;
        weight0.v = _mm256_mul_ps(weight, new_rbint0);
        weight1.v = _mm256_mul_ps(weight, rbint);

// t=0, x=1, y=0
// t=1, x=1, y=0
#define sift_calc_keypoint_descriptor_f_avx_hist_2_3(p)                        \
  if (new_binx1_biny_flag_f.f32[p]) {                                          \
    atd(new_binx1.i32[p], biny.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                        \
    atd(new_binx1.i32[p], biny.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                        \
  }
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(0)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(1)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(2)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(3)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(4)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(5)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(6)
        sift_calc_keypoint_descriptor_f_avx_hist_2_3(7)
#undef sift_calc_keypoint_descriptor_f_avx_hist_2_3
      }

      if (_mm256_movemask_ps(binx_new_biny1_flag_f.v)) {
        __m256 weight =
            _mm256_mul_ps(_mm256_mul_ps(win_mod, new_rbinx0), rbiny);
        UV256F32 weight0, weight1;
        weight0.v = _mm256_mul_ps(weight, new_rbint0);
        weight1.v = _mm256_mul_ps(weight, rbint);

// t=0, x=0, y=1
// t=1, x=0, y=1
#define sift_calc_keypoint_descriptor_f_avx_hist_4_5(p)                        \
  if (binx_new_biny1_flag_f.f32[p]) {                                          \
    atd(binx.i32[p], new_biny1.i32[p], bint.i32[p] % ORIENTATION_BINS) +=      \
        weight0.f32[p];                                                        \
    atd(binx.i32[p], new_biny1.i32[p], new_bint1.i32[p] % ORIENTATION_BINS) += \
        weight1.f32[p];                                                        \
  }
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(0)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(1)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(2)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(3)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(4)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(5)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(6)
        sift_calc_keypoint_descriptor_f_avx_hist_4_5(7)
#undef sift_calc_keypoint_descriptor_f_avx_hist_4_5
      }

      if (_mm256_movemask_ps(new_binx1_new_biny1_flag_f.v)) {
        __m256 weight = _mm256_mul_ps(_mm256_mul_ps(win_mod, rbinx), rbiny);
        UV256F32 weight0, weight1;
        weight0.v = _mm256_mul_ps(weight, new_rbint0);
        weight1.v = _mm256_mul_ps(weight, rbint);

// t=0, x=1, y=1
// t=1, x=1, y=1
#define sift_calc_keypoint_descriptor_f_avx_hist_6_7(p)                        \
  if (new_binx1_new_biny1_flag_f.f32[p]) {                                     \
    atd(new_binx1.i32[p], new_biny1.i32[p], bint.i32[p] % ORIENTATION_BINS) += \
        weight0.f32[p];                                                        \
    atd(new_binx1.i32[p], new_biny1.i32[p],                                    \
        new_bint1.i32[p] % ORIENTATION_BINS) += weight1.f32[p];                \
  }
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(0)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(1)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(2)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(3)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(4)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(5)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(6)
        sift_calc_keypoint_descriptor_f_avx_hist_6_7(7)
#undef sift_calc_keypoint_descriptor_f_avx_hist_6_7
      }

    }  // for (dxi=lb; dxi<=hb-7; dxi+=8) // avx 256

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
    }  // for (; dxi <= hb; dxi++) // c

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
