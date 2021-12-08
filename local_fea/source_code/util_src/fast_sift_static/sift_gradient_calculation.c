/*
 * sift_gradient_calculation.c
 *
 *  Created on: Nov 3, 2014
 *      Author: peter
 */

#include "sift_gradient_calculation.h"
#include "sift_math.h"
#include "sift_sse4_2_typedef.h"

void sift_update_gradient_f(SIFT_FILTER *f) {
  int x, y, s;
  const int xo = 1;
  const int yo = f->octave_width;
  const int so = f->octave_size;

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
    ;
    // middle pixels
    x = yo - 2;
    while (x--) {
      gx = 0.5 * (src[+xo] - src[-xo]);
      gy = src[+yo] - src[0];
      SAVE_BACK
      ;
    }
    // last pixel
    gx = src[0] - src[-xo];
    gy = src[+yo] - src[0];
    SAVE_BACK
    ;

    // middle rows
    for (y = 1; y < f->octave_height - 1; y++) {
      // 1st pixel
      gx = src[+xo] - src[0];
      gy = 0.5 * (src[+yo] - src[-yo]);
      SAVE_BACK
      ;
      // middle pixels
      x = yo - 2;
      while (x--) {
        gx = 0.5 * (src[+xo] - src[-xo]);
        gy = 0.5 * (src[+yo] - src[-yo]);
        SAVE_BACK
        ;
      }
      // last pixel
      gx = src[0] - src[-xo];
      gy = 0.5 * (src[+yo] - src[-yo]);
      SAVE_BACK
      ;
    }

    // last row
    // 1st pixel
    gx = src[+xo] - src[0];
    gy = src[0] - src[-yo];
    SAVE_BACK
    ;
    // middle pixels
    x = yo - 2;
    while (x--) {
      gx = 0.5 * (src[+xo] - src[-xo]);
      gy = src[0] - src[-yo];
      SAVE_BACK
      ;
    }
    // last pixel
    gx = src[0] - src[-xo];
    gy = src[0] - src[-yo];
    SAVE_BACK
    ;
  }  // for (s=f->s_min+1; s<=f->s_max-2; s++)
  f->grad_o = f->o_cur;

#undef SAVE_BACK
}
