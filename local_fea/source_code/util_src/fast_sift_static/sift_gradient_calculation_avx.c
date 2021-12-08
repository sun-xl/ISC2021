/*
 * sift_gradient_calculation_avx.c
 *
 *  Created on: Nov 10, 2014
 *      Author: peter
 */

#include "sift_gradient_calculation.h"
#include "sift_math.h"
#include "sift_avx_typedef.h"

static inline __m256 fast_sift_fast_sqrt_f_avx(const __m256 x) {
  // 32-bits version
  union {
    __m256 x;
    __m256i i;
  } u;

  const __m256 xhalf = _mm256_mul_ps(_mm256_set1_ps(0.5), x);
  const __m128i magic = _mm_set1_epi32(0x5f3759df);
  const __m256 one_half = _mm256_set1_ps(1.5);

  // convert floating point value in RAW integer
  u.x = x;

  // gives initial guess y0
  u.i = _mm256_insertf128_si256(
      _mm256_castsi128_si256(
          _mm_sub_epi32(magic, _mm_srli_epi32(_mm256_castsi256_si128(u.i), 1))),
      _mm_sub_epi32(magic, _mm_srli_epi32(_mm256_extractf128_si256(u.i, 1), 1)),
      1);

  // two Newton steps
  u.x = _mm256_mul_ps(
      u.x,
      _mm256_sub_ps(one_half, _mm256_mul_ps(xhalf, _mm256_mul_ps(u.x, u.x))));
  u.x = _mm256_mul_ps(
      u.x,
      _mm256_sub_ps(one_half, _mm256_mul_ps(xhalf, _mm256_mul_ps(u.x, u.x))));

  // zero case
  u.x = _mm256_mul_ps(
      _mm256_and_ps(u.x, _mm256_cmp_ps(_mm256_set1_ps(1e-8), x, 2)), x);
  return u.x;
}

// use _mm256_atan_ps
static inline __m256 fast_atan2_mode_mod_2pi_f_avx(const __m256 y256,
                                                   const __m256 x256) {
  const __m256 PI2v = _mm256_set1_ps((float)PI2);
  const __m256 ZEROv = _mm256_set1_ps(0.f);
  __m256 absy =
      _mm256_add_ps(_mm256_andnot_ps(_mm256_set1_ps(-0.f), y256),  // fabs
                    _mm256_set1_ps(EPSILON_F));
  __m256 xflag = _mm256_cmp_ps(x256, ZEROv, 13);  // if x>=0, 1s; else 0s
  __m256 yflag = _mm256_cmp_ps(y256, ZEROv, 13);  // if y>=0, 1s; else 0s
  __m256 x_plus_absy = _mm256_add_ps(x256, absy);

  __m256 r = _mm256_div_ps(
      _mm256_or_ps(_mm256_and_ps(xflag, _mm256_sub_ps(x256, absy)),
                   _mm256_andnot_ps(xflag, x_plus_absy)),
      _mm256_or_ps(_mm256_and_ps(xflag, x_plus_absy),
                   _mm256_andnot_ps(xflag, _mm256_sub_ps(absy, x256))));

  __m256 angle = _mm256_or_ps(
      _mm256_and_ps(xflag, _mm256_set1_ps((float)(PI / 4))),
      _mm256_andnot_ps(xflag, _mm256_set1_ps((float)(3 * PI / 4))));

  angle = _mm256_add_ps(
      angle, _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(0.1821F),
                                                       _mm256_mul_ps(r, r)),
                                         _mm256_set1_ps(0.9675F)),
                           r));

  angle = _mm256_or_ps(_mm256_and_ps(yflag, angle),
                       _mm256_andnot_ps(yflag, _mm256_sub_ps(ZEROv, angle)));

  angle = _mm256_add_ps(PI2v, angle);

  while (1) {
    __m256 cmp_res = _mm256_cmp_ps(angle, PI2v, 14);  // gt
    int mask = _mm256_movemask_ps(cmp_res);
    if (!mask) break;
    angle = _mm256_sub_ps(angle, _mm256_and_ps(cmp_res, PI2v));
  }
  while (1) {
    __m256 cmp_res = _mm256_cmp_ps(angle, ZEROv, 1);  // lt
    int mask = _mm256_movemask_ps(cmp_res);
    if (!mask) break;
    angle = _mm256_add_ps(angle, _mm256_and_ps(cmp_res, PI2v));
  }

  return angle;
}

void sift_update_gradient_f_avx(SIFT_FILTER *f) {
  int x, y, s;
  const int xo = 1;
  const int yo = f->octave_width;
  const int so = f->octave_size;
  const __m256 half = _mm256_set1_ps(0.5);

  if (f->grad_o == f->o_cur)  // no need to update gradient
    return;

  // GSS planar [s_min, s_max]
  // gradient planar, DoG planar [s_min+1, s_max-2]
  for (s = f->s_min + 1; s <= f->s_max - 2; s++)  // for each planar
  {
    float *src, *grad, gx, gy;

// save magnitude & orientation
#define SAVE_BACK                                               \
  *grad++ = sift_fast_sqrt_f(gx * gx + gy * gy);                \
  *grad++ = sift_mod_2pi_f(sift_fast_atan2_f(gy, gx) + 2 * PI); \
  ++src;

    src = f->octave + (s - f->s_min) * so;
    grad = f->grad + (s - f->s_min - 1) * so * 2;  // stride = 2*so

    // 1st row
    // 1st pixel
    gx = src[+xo] - src[0];
    gy = src[+yo] - src[0];
    SAVE_BACK
    // middle pixels
    x = yo - 2;
    while (x--) {
      gx = 0.5 * (src[+xo] - src[-xo]);
      gy = src[+yo] - src[0];
      SAVE_BACK
    }
    // last pixel
    gx = src[0] - src[-xo];
    gy = src[+yo] - src[0];
    SAVE_BACK

    // middle rows
    for (y = 1; y < f->octave_height - 1; y++) {
      // 1st pixel
      gx = src[+xo] - src[0];
      gy = 0.5 * (src[+yo] - src[-yo]);
      SAVE_BACK

      // middle pixels
      const int flagxp = ((uintptr_t)(src - xo)) & 0x1F;
      ;
      const int flagxn = ((uintptr_t)(src + xo)) & 0x1F;
      const int flagyp = ((uintptr_t)(src - yo)) & 0x1F;
      const int flagyn = ((uintptr_t)(src + yo)) & 0x1F;
      const int flag = ((uintptr_t)(grad)) & 0x1F;
      __m256 srcp256, srcn256, gx256, gy256;
      for (x = yo - 2; x >= 8; x -= 8) {
        if (flagxp)
          srcp256 = _mm256_loadu_ps(src - xo);
        else
          srcp256 = _mm256_load_ps(src - xo);

        if (flagxn)
          srcn256 = _mm256_loadu_ps(src + xo);
        else
          srcn256 = _mm256_load_ps(src + xo);
        gx256 = _mm256_mul_ps(half, _mm256_sub_ps(srcn256, srcp256));

        if (flagyp)
          srcp256 = _mm256_loadu_ps(src - yo);
        else
          srcp256 = _mm256_load_ps(src - yo);
        if (flagyn)
          srcn256 = _mm256_loadu_ps(src + yo);
        else
          srcn256 = _mm256_load_ps(src + yo);
        gy256 = _mm256_mul_ps(half, _mm256_sub_ps(srcn256, srcp256));

        // magnitude & orientation
        srcp256 = fast_sift_fast_sqrt_f_avx(_mm256_add_ps(
            _mm256_mul_ps(gx256, gx256), _mm256_mul_ps(gy256, gy256)));
        srcn256 = fast_atan2_mode_mod_2pi_f_avx(gy256, gx256);

        // srcp256 = _mm256_sqrt_ps (_mm256_add_ps (_mm256_mul_ps (gx256,
        // gx256), _mm256_mul_ps (gy256, gy256)));
        // srcn256 = _mm256_atan2_ps (gy256, gx256);

        // m7-m6-m5-m4-m3-m2-m1-m0
        // o7-o6-o5-o4-o3-o2-o1-o0
        gx256 =
            _mm256_unpacklo_ps(srcp256, srcn256);  // o5-m5-o4-m4-o1-m1-o0-m0
        gy256 =
            _mm256_unpackhi_ps(srcp256, srcn256);  // o7-m7-o6-m6-o3-m3-o2-m2

        srcp256 = _mm256_insertf128_ps(gx256, _mm256_castps256_ps128(gy256),
                                       1);  // o3-m3-o2-m2-o1-m1-o0-m0
        srcn256 = _mm256_insertf128_ps(gy256, _mm256_extractf128_ps(gx256, 1),
                                       0);  // o7-m7-o6-m6-o5-m5-o4-m4

        if (flag) {
          _mm256_storeu_ps(grad, srcp256);
          grad += 8;
          _mm256_storeu_ps(grad, srcn256);
          grad += 8;
        } else {
          _mm256_store_ps(grad, srcp256);
          grad += 8;
          _mm256_store_ps(grad, srcn256);
          grad += 8;
        }
        src += 8;
      }

      while (x--) {
        gx = 0.5 * (src[+xo] - src[-xo]);
        gy = 0.5 * (src[+yo] - src[-yo]);
        SAVE_BACK
      }
      // last pixel
      gx = src[0] - src[-xo];
      gy = 0.5 * (src[+yo] - src[-yo]);
      SAVE_BACK
    }

    // last row
    // 1st pixel
    gx = src[+xo] - src[0];
    gy = src[0] - src[-yo];
    SAVE_BACK
    // middle pixels
    x = yo - 2;
    while (x--) {
      gx = 0.5 * (src[+xo] - src[-xo]);
      gy = src[0] - src[-yo];
      SAVE_BACK
    }
    // last pixel
    gx = src[0] - src[-xo];
    gy = src[0] - src[-yo];
    SAVE_BACK
  }  // for (s=f->s_min+1; s<=f->s_max-2; s++)
  f->grad_o = f->o_cur;

#undef SAVE_BACK
}
