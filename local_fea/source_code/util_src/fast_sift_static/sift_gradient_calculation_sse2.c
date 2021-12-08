/*
 * sift_gradient_calculation_sse2.c
 *
 *  Created on: Nov 10, 2014
 *      Author: peter
 */

#include "sift_gradient_calculation.h"
#include "sift_math.h"
#include "sift_sse4_2_typedef.h"

static inline __m128 fast_sift_fast_sqrt_f_sse2(const __m128 x) {
  // 32-bit version
  union {
    __m128 f32_v128;
    __m128i i32_v128;
  } u;
  // constant
  const __m128 one_half_v128 = _mm_set1_ps(1.5f);
  const __m128 xhalf = _mm_mul_ps(_mm_set1_ps(0.5f), x);

  // convert floating point value in RAW integer
  u.f32_v128 = x;

  // gives initial guess y0
  u.i32_v128 =
      _mm_sub_epi32(_mm_set1_epi32(0x5f3759df), _mm_srli_epi32(u.i32_v128, 1));

  // two Newton steps
  u.f32_v128 = _mm_mul_ps(
      u.f32_v128,
      _mm_sub_ps(one_half_v128,
                 _mm_mul_ps(xhalf, _mm_mul_ps(u.f32_v128, u.f32_v128))));
  u.f32_v128 = _mm_mul_ps(
      u.f32_v128,
      _mm_sub_ps(one_half_v128,
                 _mm_mul_ps(xhalf, _mm_mul_ps(u.f32_v128, u.f32_v128))));

  // zero case
  u.f32_v128 =
      _mm_mul_ps(x, _mm_and_ps(u.f32_v128, _mm_cmpge_ps(x, _mm_set1_ps(1e-8))));
  return u.f32_v128;
}

static inline __m128 fast_atan2_mode_mod_2pi_f_sse2(const __m128 y,
                                                    const __m128 x) {
  const __m128 ZERO_V128 = _mm_set1_ps(0.f);
  const __m128 PI2_V128 = _mm_set1_ps(PI2);

  const __m128 abs_y =
      _mm_add_ps(_mm_andnot_ps(_mm_set1_ps(-0.f), y), _mm_set1_ps(EPSILON_F));
  const __m128 x_abs_y = _mm_add_ps(x, abs_y);
  const __m128 x_flag = _mm_cmpge_ps(x, ZERO_V128);
  const __m128 y_flag = _mm_cmpge_ps(y, ZERO_V128);

  __m128 r = _mm_div_ps(_mm_or_ps(_mm_and_ps(x_flag, _mm_sub_ps(x, abs_y)),
                                  _mm_andnot_ps(x_flag, x_abs_y)),
                        _mm_or_ps(_mm_and_ps(x_flag, x_abs_y),
                                  _mm_andnot_ps(x_flag, _mm_sub_ps(abs_y, x))));

  __m128 angle;
  angle = _mm_or_ps(_mm_and_ps(x_flag, _mm_set1_ps((float)PI / 4)),
                    _mm_andnot_ps(x_flag, _mm_set1_ps((float)(3 * PI / 4))));

  angle = _mm_add_ps(
      angle,
      _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(_mm_set1_ps(0.1821F), _mm_mul_ps(r, r)),
                            _mm_set1_ps(0.9675F)),
                 r));

  angle = _mm_or_ps(_mm_and_ps(y_flag, angle),
                    _mm_andnot_ps(y_flag, _mm_sub_ps(ZERO_V128, angle)));

  angle = _mm_add_ps(PI2_V128, angle);

  // clamping
  while (1) {
    __m128 cmp_res = _mm_cmpgt_ps(angle, PI2_V128);
    int mask = _mm_movemask_ps(cmp_res);
    if (!mask) break;
    angle = _mm_sub_ps(angle, _mm_and_ps(cmp_res, PI2_V128));
  }
  while (1) {
    __m128 cmp_res = _mm_cmplt_ps(angle, ZERO_V128);
    int mask = _mm_movemask_ps(cmp_res);
    if (!mask) break;
    angle = _mm_add_ps(angle, _mm_and_ps(cmp_res, PI2_V128));
  }

  return angle;
}

void sift_update_gradient_f_sse2(SIFT_FILTER *f) {
  int x, y, s;
  const int xo = 1;
  const int yo = f->octave_width;
  const int so = f->octave_size;
  const __m128 half_f32 = _mm_set1_ps(0.5f);

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
      // 1: unaligned, 0: aligned
      const unsigned int prev_x_align_flag = ((uintptr_t)(src - xo)) & 0x0F;
      const unsigned int next_x_align_flag = ((uintptr_t)(src + xo)) & 0x0F;
      const unsigned int prev_y_align_flag = ((uintptr_t)(src - yo)) & 0x0F;
      const unsigned int next_y_align_flag = ((uintptr_t)(src + yo)) & 0x0F;
      const unsigned int grad_align_flag = ((uintptr_t)(grad)) & 0xF;
      __m128 src_prev_v128, src_next_v128, gradx_v128, grady_v128;
      __m128 magnitude_v128, orientation_v128;
      for (x = yo - 2; x >= 4; x -= 4)  // sse 128
      {
        if (prev_x_align_flag)  // unaligned
          src_prev_v128 = _mm_loadu_ps(src - xo);
        else  // aligned
          src_prev_v128 = _mm_load_ps(src - xo);

        if (next_x_align_flag)  // unaligned
          src_next_v128 = _mm_loadu_ps(src + xo);
        else
          src_next_v128 = _mm_load_ps(src + xo);

        gradx_v128 =
            _mm_mul_ps(half_f32, _mm_sub_ps(src_next_v128, src_prev_v128));

        if (prev_y_align_flag)  // unaligned
          src_prev_v128 = _mm_loadu_ps(src - yo);
        else  // aligned
          src_prev_v128 = _mm_load_ps(src - yo);

        if (next_y_align_flag)  // unaligned
          src_next_v128 = _mm_loadu_ps(src + yo);
        else  // aligned
          src_next_v128 = _mm_load_ps(src + yo);

        grady_v128 =
            _mm_mul_ps(half_f32, _mm_sub_ps(src_next_v128, src_prev_v128));

        // magnitude & orientation
        // m3-m2-m1-m0
        // o3-o2-o1-o0
        magnitude_v128 = fast_sift_fast_sqrt_f_sse2(
            _mm_add_ps(_mm_mul_ps(gradx_v128, gradx_v128),
                       _mm_mul_ps(grady_v128, grady_v128)));
        orientation_v128 =
            fast_atan2_mode_mod_2pi_f_sse2(grady_v128, gradx_v128);

        if (grad_align_flag)  // unaligned
        {
          _mm_storeu_ps(grad,
                        _mm_unpacklo_ps(magnitude_v128, orientation_v128));
          grad += 4;
          _mm_storeu_ps(grad,
                        _mm_unpackhi_ps(magnitude_v128, orientation_v128));
          grad += 4;
        } else  // aligned
        {
          _mm_store_ps(grad, _mm_unpacklo_ps(magnitude_v128, orientation_v128));
          grad += 4;
          _mm_store_ps(grad, _mm_unpackhi_ps(magnitude_v128, orientation_v128));
          grad += 4;
        }

        src += 4;
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
