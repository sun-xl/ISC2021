/*
 * sift_keypoint_detection.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_keypoint_detection.h"
#include "sift_sse4_2_typedef.h"
#include <stdlib.h>
#include <math.h>

static void gaussian_elimination_1x3_3x3_f(float A[3 * 3], float b[3]) {
  int i, j;
  int row, col;
  float tmp;

  for (j = 0; j < 3; j++) {
    float maxa = 0.;
    float maxabsa = 0.;
    int maxi = -1;

    // look for the maximally stable pivot
    for (i = j; i < 3; i++) {
      float a = A[j * 3 + i];
      float absa = fabsf(a);
      if (absa > maxabsa) {
        maxa = a;
        maxabsa = absa;
        maxi = i;
      }
    }  // for (i=j; i<3; i++)

    // if singular give up
    if (maxabsa < 1e-10f) {
      b[0] = 0;
      b[1] = 0;
      b[2] = 0;
      break;
    }

    // swap j-the column with maxi-th column, then normalize j-th column
    for (row = j; row < 3; row++) {
      tmp = A[row * 3 + maxi];
      A[row * 3 + maxi] = A[row * 3 + j];
      A[row * 3 + j] = tmp;
      A[row * 3 + j] /= maxa;
    }
    tmp = b[maxi];
    b[maxi] = b[j];
    b[j] = tmp;
    b[j] /= maxa;

    // use j-th column to eliminate the following columns
    for (col = j + 1; col < 3; col++)  // the following columns
    {
      tmp = A[j * 3 + col];
      for (row = j; row < 3; row++) A[row * 3 + col] -= tmp * A[row * 3 + j];
      b[col] -= tmp * b[j];
    }
  }  // for (j=0; j<3; j++)

  // 							 |1,   0,   0  |
  // |deltax, deltay, deltas|. |??,  1,   0  | = |??, ??, ??|
  // 							 |??,  ??,  1  |

  // 							 |1,   0,   0  |
  // |deltax, deltay, deltas|. |0 ,  1,   0  | = |??, ??, ??|
  // 							 |0 ,  0 ,  1  |
  // backward substitution
  for (i = 2; i > 0; i--) {
    tmp = b[i];
    for (col = i - 1; col >= 0; col--)  // the frontal columns
      b[col] -= tmp * A[i * 3 + col];   // j==i
  }
}

void sift_detect_f(SIFT_FILTER *f) {
  int x, y, s, k;
  float *dog;
  float v;
  SIFT_KEYPOINT *kc;
  const int xo = 1;
  const int yo = f->octave_width;
  const int so = f->octave_size;

  f->nkeys = 0;  // important!!!

  // compute DoG
  dog = f->dog;
  for (s = f->s_min; s <= f->s_max - 1; s++) {
    int size = so;
    float *sc, *sn;
    sc = f->octave + (s - f->s_min) * so;
    sn = sc + so;
    while (size--) *dog++ = *sn++ - *sc++;
  }  // for (s=f->s_min; s<=f->s_max-1; s++)

  // find local extrema of DoG
  // considering boundary, starting from DoG [1,1,f->s_min+1]
  dog = f->dog + xo + yo + so;
  for (s = f->s_min + 1; s <= f->s_max - 2; s++) {
    for (y = 1; y < f->octave_height - 1; y++) {
      for (x = 1; x < f->octave_width - 1; x++) {
        v = *dog;
// local extrama detector & low contrast threshold
#define CHECK_NEIGHBORS(CMP, SGN)                                       \
  (v CMP## = SGN 0.8 * f->peak_thresh && v CMP * (dog - yo - xo) &&     \
             v CMP * (dog - yo) && v CMP * (dog - yo + xo) &&           \
             v CMP * (dog - xo) && v CMP * (dog + xo) &&                \
             v CMP * (dog + yo - xo) && v CMP * (dog + yo) &&           \
             v CMP * (dog + yo + xo) && v CMP * (dog - yo - xo + so) && \
             v CMP * (dog - yo + so) && v CMP * (dog - yo + xo + so) && \
             v CMP * (dog - xo + so) && v CMP * (dog + so) &&           \
             v CMP * (dog + xo + so) && v CMP * (dog + yo - xo + so) && \
             v CMP * (dog + yo + so) && v CMP * (dog + yo + xo + so) && \
             v CMP * (dog - yo - xo - so) && v CMP * (dog - yo - so) && \
             v CMP * (dog - yo + xo - so) && v CMP * (dog - xo - so) && \
             v CMP * (dog - so) && v CMP * (dog + xo - so) &&           \
             v CMP * (dog + yo - xo - so) && v CMP * (dog + yo - so) && \
             v CMP * (dog + yo + xo - so))
        if (CHECK_NEIGHBORS(>, +) || CHECK_NEIGHBORS(<, -)) {
          // make space for more keypoints
          if (f->nkeys >= f->keys_res) {
            f->keys_res += 500;
            f->keys = realloc(f->keys, f->keys_res * sizeof(SIFT_KEYPOINT));
          }
          kc = f->keys + (f->nkeys++);
          kc->ix = x;
          kc->iy = y;
          kc->is = s;
        }           // if (CHECK_NEIGHBORS(>,+) || CHECK_NEIGHBORS(<,-))
        dog++;      // next pixel
      }             // for (x=1; x<f->octave_width-1; x++)
      dog += 2;     // next row
    }               // for (y=1; y<f->octave_height; y++)
    dog += 2 * yo;  // next planar
  }                 // for (s=f->s_min+1; s<=f->s_max-2; s++)

  // accurate keypoint localization - refine the local extrama
  kc = f->keys;
  for (k = 0; k < f->nkeys; k++)  // for each keypoint
  {
    float Dx, Dy, Ds;
    float Dxx, Dyy, Dss;
    float Dxy, Dxs, Dys;
    float A[3 * 3], b[3];
    int dx, dy, iter;

    dx = 0;
    dy = 0;
    x = f->keys[k].ix;
    y = f->keys[k].iy;
    s = f->keys[k].is;
    for (iter = 0; iter < 5; iter++) {
      // prepare for next iteration
      x += dx;
      y += dy;

      dog = f->dog + x * xo + y * yo + (s - f->s_min) * so;
#define at(dx, dy, ds) \
  (*(dog + (dx) * xo + (dy) * yo + (ds) * so))  // DoG manipulation
      // compute the gradient
      Dx = 0.5 * (at(+1, 0, 0) - at(-1, 0, 0));
      Dy = 0.5 * (at(0, +1, 0) - at(0, -1, 0));
      Ds = 0.5 * (at(0, 0, +1) - at(0, 0, -1));

      // compute the Hessian
      Dxx = (at(+1, 0, 0) + at(-1, 0, 0) - 2.0 * at(0, 0, 0));
      Dyy = (at(0, +1, 0) + at(0, -1, 0) - 2.0 * at(0, 0, 0));
      Dss = (at(0, 0, +1) + at(0, 0, -1) - 2.0 * at(0, 0, 0));

      Dxy = 0.25 *
            (at(+1, +1, 0) + at(-1, -1, 0) - at(-1, +1, 0) - at(+1, -1, 0));
      Dxs = 0.25 *
            (at(+1, 0, +1) + at(-1, 0, -1) - at(-1, 0, +1) - at(+1, 0, -1));
      Dys = 0.25 *
            (at(0, +1, +1) + at(0, -1, -1) - at(0, -1, +1) - at(0, +1, -1));

      // solve linear system
      // 							 |Dxx, Dxy, Dxs|
      // |deltax, deltay, deltas|. |Dxy, Dyy, Dys| = |-Dx, -Dy, -Ds|
      // 							 |Dxs, Dys, Dss|
      A[0] = Dxx;
      A[3 * 1 + 1] = Dyy;
      A[3 * 2 + 2] = Dss;
      A[3 * 1] = A[1] = Dxy;
      A[3 * 2] = A[2] = Dxs;
      A[3 * 2 + 1] = A[3 * 1 + 2] = Dys;

      b[0] = -Dx;
      b[1] = -Dy;
      b[2] = -Ds;

      // Gaussian elimination
      gaussian_elimination_1x3_3x3_f(A, b);

      // If the translation of the keypoint is big, move the keypoint
      // and re-iterate the computation. Otherwise we are all set.
      dx = ((b[0] > 0.6 && x < f->octave_width - 2) ? 1 : 0) +
           ((b[0] < -0.6 && x > 1) ? -1 : 0);

      dy = ((b[1] > 0.6 && y < f->octave_height - 2) ? 1 : 0) +
           ((b[1] < -0.6 && y > 1) ? -1 : 0);
      if (dx == 0 && dy == 0) break;
    }  // for (iter=0; iter<5; iter++)

    // check threshold and other conditions
    {
      float val = at(0, 0, 0) + 0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2]);
      float score = (Dxx + Dyy) * (Dxx + Dyy) / (Dxx * Dyy - Dxy * Dxy);
      float xn = x + b[0];
      float yn = y + b[1];
      float sn = s + b[2];

      int good = fabsf(val) > f->peak_thresh &&
                 score < (f->edge_thresh + 1) * (f->edge_thresh + 1) /
                             f->edge_thresh &&
                 score >= 0 && fabsf(b[0]) < 1.5 && fabsf(b[1]) < 1.5 &&
                 fabsf(b[2]) < 1.5 && xn >= 0 && xn <= f->octave_width - 1 &&
                 yn >= 0 && yn <= f->octave_height - 1 && sn >= f->s_min &&
                 sn <= f->s_max;

      if (good) {
        kc->o = f->o_cur;
        kc->ix = x;
        kc->iy = y;
        kc->is = s;
        kc->s = sn;
        kc->x = xn * f->xper;
        kc->y = yn * f->xper;
        kc->sigma = f->sigma0 * pow(2.0, sn / f->S) * f->xper;
        ++kc;
      }
    }  // done checking threshold and other conditions
  }    // for (k=0; k<f->nkeys; k++)
       // update keypoint count
  f->nkeys = (int)(kc - f->keys);

#undef CHECK_NEIGHBORS
#undef at
}
