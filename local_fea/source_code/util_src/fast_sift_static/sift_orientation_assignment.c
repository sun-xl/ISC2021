/*
 * sift_orientation_assignment.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_orientation_assignment.h"
#include "sift_math.h"

#include <string.h>

#define SIFT_BILINEAR_ORIENTATIONS

int sift_calc_keypoint_orientations_f(const SIFT_FILTER *const f,
                                      float angles[4],
                                      const SIFT_KEYPOINT *const k) {
  const int xo = 2;  // magnitude & orientation
  const int yo = 2 * f->octave_width;
  const int so = 2 * f->octave_size;

  const float x = k->x / f->xper;
  const float y = k->y / f->xper;
  const float sigma = k->sigma / f->xper;

  const int xi = (int)(x + 0.5);
  const int yi = (int)(y + 0.5);
  const int si = k->is;

  const float *pt = f->grad + xo * xi + yo * yi + so * (si - f->s_min - 1);

  const float winf = 1.5f;
  const float sigmaW = winf * sigma;
  const int W = SIFT_MAX(floorf(3.0 * sigmaW), 1);  // Gaussian window width

  const int nbins = 36;  // 2*pi/36
  float hist[nbins];

  int xs, ys;
  int iter, i;
  float maxh;
  int nangles;

  if (k->o != f->o_cur)  // skip if the keypoint octave is not current
    return 0;

  // skip the keypoint if it is out of bounds
  if (xi < 0 || xi > f->octave_width - 1 || yi < 0 ||
      yi > f->octave_height - 1 || si < f->s_min + 1 || si > f->s_max - 2)
    return 0;

  // clear histogram
  memset(hist, 0, sizeof(float) * nbins);

  // fill histogram
  for (ys = SIFT_MAX(-W, -yi); ys <= SIFT_MIN(W, f->octave_height - 1 - yi);
       ys++) {
    for (xs = SIFT_MAX(-W, -xi); xs <= SIFT_MIN(W, f->octave_width - 1 - xi);
         xs++) {
      float dx = (float)(xi + xs) - x;
      float dy = (float)(yi + ys) - y;
      float r2 = dx * dx + dy * dy;
      float wgt, mod, ang, fbin;

      // limit to a circular window
      if (r2 >= W * W + 0.6) continue;

      wgt = fast_expn_f(r2 / (2 * sigmaW * sigmaW));  // Gaussian weight
      mod = *(pt + xs * xo + ys * yo);                // magnitude
      ang = *(pt + xs * xo + ys * yo + 1);            // orientation
      fbin = nbins * ang / (2 * PI);                  // float bin index

#if defined(SIFT_BILINEAR_ORIENTATIONS)
      {
        int bin = (int)sift_floor_f(fbin - 0.5);
        float rbin = fbin - bin - 0.5f;
        hist[(bin + nbins) % nbins] += (1 - rbin) * mod * wgt;
        hist[(bin + 1) % nbins] += (rbin) * mod * wgt;
      }
#else
      {
        int bin = sift_floor_d(fbin);
        hist[(bin) % nbins] += mod * wgt;
      }
#endif
    }  // for (ys=MAX(-W, -yi); ys<=MIN(W, f->octave_height-1-yi); ys++)
  }    // for (ys=MAX(-W, -yi); ys<=MIN(W, f->octave_height-1-yi); ys++)

  // smooth histogram
  for (iter = 0; iter < 6; iter++) {
    float prev = hist[nbins - 1];
    float first = hist[0];
    int i;
    for (i = 0; i < nbins - 1; i++) {
      float newh = (prev + hist[i] + hist[(i + 1) % nbins]) / 3.0;  // mean
      prev = hist[i];
      hist[i] = newh;
    }
    hist[i] = (prev + hist[i] + first) / 3.0;  // i == nbins - 1
  }

  // find the histogram maximum
  maxh = 0;
  for (i = 0; i < nbins; ++i) maxh = SIFT_MAX(maxh, hist[i]);

  // find peaks within 80% from max
  nangles = 0;
  for (i = 0; i < nbins; ++i) {
    float h0 = hist[i];
    float hm = hist[(i - 1 + nbins) % nbins];
    float hp = hist[(i + 1 + nbins) % nbins];

    // is this a peak
    if (h0 > 0.8 * maxh && h0 > hm && h0 > hp) {
      // quadratic interpolation
      float di = -0.5 * (hp - hm) / (hp + hm - 2 * h0);
      float th = 2 * PI * (i + di + 0.5) / nbins;
      angles[nangles++] = th;
      if (nangles == 4) break;
    }
  }

  return nangles;
}
